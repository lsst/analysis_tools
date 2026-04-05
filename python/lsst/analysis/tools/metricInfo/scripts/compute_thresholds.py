#!/usr/bin/env python3
# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Compute MAD-based per-stratum thresholds from a drp_pipe pipeline.

Parses a drp_pipe pipeline YAML to discover all metric-producing tasks and
their corresponding metric tables (produced by MakeMetricTableTask).  Checks
that every metric bundle type is covered by a MakeMetricTableTask, and that
every expected table dataset actually exists in the butler collection.
Thresholds are computed only for tables that are present; missing tables are
reported so the user can generate them before re-running.

Outputs:
- A threshold CSV mirroring the pipeline path structure:
    $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml → thresholds/LSSTCam/DRP.csv
- A PDF report showing histograms of each metric's distribution per stratum,
  with threshold lines overlaid, grouped by metricType:
    thresholds/LSSTCam/DRP.pdf

Examples
--------
::

    python compute_thresholds.py \\
        --butler-repo /repo/main \\
        --collection LSSTCam/runs/DRP/20250501_20250609/w_2025_30/DM-51933 \\
        --pipeline $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml
"""
from __future__ import annotations

import argparse
import logging
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from generate_report import generate_report
from astropy.coordinates import SkyCoord
import astropy.units as u

from lsst.daf.butler import Butler, Config
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.pipe.base import Pipeline

log = logging.getLogger(__name__)

_METRICINFO_DIR = Path(__file__).resolve().parent.parent
METRIC_DESCRIPTIONS_PATH = _METRICINFO_DIR / "metricDescriptions.yaml"
THRESHOLDS_DIR = _METRICINFO_DIR / "thresholds"


# MAD_SCALE is a consistency factor that makes MAD comparable to a standard
# deviation for normally distributed data (MAD_SCALE * MAD ≈ σ).  Keeps
# --n-mad intuitive as an "n-sigma equivalent" threshold width.
MAD_SCALE = 1.4826
DEFAULT_N_MAD = 5.0
DEFAULT_TIER1_FRACTION = 0.05
MIN_VISITS_PER_STRATUM = 10


def load_metric_properties(path: Path) -> dict[tuple[str, str, str], dict]:
    """Load metricDescriptions.yaml and flatten to a lookup dict.

    Parameters
    ----------
    path
        Path to metricDescriptions.yaml.

    Returns
    -------
    dict
        Keyed by ``(task_class, atool, metric_name)``; values are the
        per-metric property dicts (``sided``, ``units``, ``metricTypes``, etc.).
    """
    with open(path) as file:
        md = yaml.safe_load(file)
    metric_dict: dict[tuple[str, str, str], dict] = {}
    for task_class, task_cfg in md["tasks"].items():
        for atool, metric in task_cfg["atools"].items():
            for metric_name, metric_properties in metric.items():
                metric_dict[(task_class, atool, metric_name)] = metric_properties
    return metric_dict


def load_bins(instrument: str) -> dict:
    """Load the instrument-specific bin definitions.

    Parameters
    ----------
    instrument
        Instrument name matching a subdirectory of ``thresholds/``.

    Returns
    -------
    dict
        Parsed bins.yaml content.
    """
    bins_path = THRESHOLDS_DIR / instrument / "bins.yaml"
    with open(bins_path) as file:
        return yaml.safe_load(file)


def extract_task_metrics(config_dict: dict) -> dict[str, dict[str, list[str]]]:
    """Extract atool, metric name template, and band list from a task config dict.

    Metric name templates differ from metric names in that they may contain ``{band}``.
    The metric template name is preserved as the key and the bands from the task's ``bands``
    config field are stored as the value.  For non-banded metrics the value is an empty list.
    This allows the metric name to be constructed from the template plus the bands.
    
    Parameters
    ----------
    config_dict
        Result of ``task_node.config.toDict()``.

    Returns
    -------
    dict
        dict of dicts. Top level dict is keyed by atool name; values are dicts keyed by metric
        template name with values a list of bands (empty list for non-banded metrics).
    """
    task_metrics: dict[str, dict[str, list[str]]] = {}
    bands: list[str] = config_dict.get("bands", [])
    for atool, atool_cfg in config_dict.get("atools", {}).items():
        produce_cfg = atool_cfg.get("produce", {})
        metric_cfg = produce_cfg.get("metric", {})
        if "units" not in metric_cfg:
            continue
        new_names = metric_cfg.get("newNames", {})
        metric_names = [new_names.get(n, n) for n in metric_cfg["units"]]
        metrics: dict[str, list[str]] = {}
        for name in metric_names:
            if "{band}" in name:
                if bands:
                    metrics[name] = list(bands)
                else:
                    # This should never happen, but just in case:
                    log.warning(
                        "Metric '%s/%s' contains {band} but task has no bands configured — skipping.",
                        atool, name,
                    )
            else:
                metrics[name] = []
        if metrics:
            task_metrics[atool] = metrics
    return task_metrics


def walk_pipeline(pipeline_path: str, butler: Butler) -> dict[str, dict]:
    """Walk through a drp_pipe pipeline and return metric table info.

    Walks the pipeline graph to find all analysis tasks producing metric
    bundles, then checks which of those bundle types are consumed by a
    MakeMetricTableTask.  Bundle types with no corresponding table task are
    reported as warnings.

    Parameters
    ----------
    pipeline_path
        Path to the pipeline YAML file.
    butler
        Butler instance used for pipeline graph resolution.

    Returns
    -------
    dict
        Keyed by metric table dataset type name.  Each value is a dict with
        keys ``task_name``, ``task_class``, ``atools`` (atool → [metric_names]),
        and ``bundle_dataset_type``.
    """
    pipeline = Pipeline.from_uri(pipeline_path)
    pipeline_graph = pipeline.to_graph(registry=butler.registry, visualization_only=True)

    # Map bundle dataset type to analysis task info
    bundle_to_task: dict[str, dict] = {}
    for task_name, task in pipeline_graph.tasks.items():
        if "metrics" not in task.outputs:
            continue
        bundle_name = task.outputs["metrics"].dataset_type_name
        atools = extract_task_metrics(task.config.toDict())
        if atools:
            bundle_to_task[bundle_name] = {
                "task_name": task_name,
                "task_class": task.task_class_name,
                "atools": atools,
            }

    # Map bundle dataset type to table dataset type via MakeMetricTableTask
    bundle_to_table: dict[str, str] = {}
    for task_name, task in pipeline_graph.tasks.items():
        if "MakeMetricTableTask" not in task.task_class_name:
            continue
        bundle_name = task.inputs["data"].dataset_type_name
        table_name = task.outputs["metricTable"].dataset_type_name
        bundle_to_table[bundle_name] = table_name

    # Warn about metric bundles that have not been fed into MakeMetricTableTask — 
    # these metrics have not been tabulated, so will not have thresholds computed.
    uncovered = set(bundle_to_task) - set(bundle_to_table)
    if uncovered:
        log.warning(
            "The following metric bundle types have no MakeMetricTableTask configured "
            "and will be skipped:\n  %s",
            "\n  ".join(sorted(uncovered)),
        )

    table_provenance: dict[str, dict] = {}
    for bundle_name, table_name in bundle_to_table.items():
        table_provenance[table_name] = {
            "task_name": bundle_to_task[bundle_name]["task_name"],
            "task_class": bundle_to_task[bundle_name]["task_class"],
            "atools": bundle_to_task[bundle_name]["atools"], # Contains metrics
            "bundle_dataset_type": bundle_name,
        }

    return table_provenance


# Compute galactic latitude for binning
def galactic_latitude(ra_deg: float, dec_deg: float) -> float:
    """Return galactic latitude in degrees for a given pointing."""
    coord = SkyCoord(ra=ra_deg * u.degree, dec=dec_deg * u.degree, frame="icrs")
    return float(coord.galactic.b.deg)


def assign_bin(value: float, edges: list[float], labels: list[str]) -> str:
    """Assign a label to a bin, which is defined by "edges".

    Parameters
    ----------
    value
        Value to bin.
    edges
        Monotonically increasing bin edges (N+1 values for N bins).
    labels
        Bin labels (N values).  Values above the last edge clamp to
        ``labels[-1]``.
    """
    for label, lo, hi in zip(labels, edges[:-1], edges[1:]):
        if lo <= value < hi:
            return label
    # If it hasn't returned yet, then it must be in the last bin.
    return labels[-1]


# Conditions lookup
def load_visit_properties(butler: Butler, collection: str) -> dict[int, dict]:
    """Build a visit property lookup dict from the visit table dataset.

    Loads the ``visitTable`` dataset once (one row per visit) and extracts
    boresight RA/Dec to compute galactic latitude.
    Further properties (e.g. PSF FWHM, moon phase) could be added here should
    more bins be needed.

    Parameters
    ----------
    butler
        Initialised Butler.
    collection
        Collection to search for the visit table dataset.

    Returns
    -------
    dict
        Keyed by visit ID; values have key ``gal_lat``.
    """
    # TO-DO: "visit_table" is specific to some pipelines. Other pipelines use a
    # different naming convention. This should either be a parameter, or extracted
    # from the pipeline.
    refs = list(butler.registry.queryDatasets("visit_table", collections=collection))
    if not refs:
        log.warning("No visit_table dataset found — galactic latitude will be 'unknown'.")
        return {}

    visit_table = butler.get(refs[0])
    log.info("Loaded visit_table with %d rows.", len(visit_table))

    # Build a lookup dict for quick referencing:
    return {
        int(row["visit"]): {"gal_lat": galactic_latitude(row["ra"], row["dec"])}
        for row in visit_table
    }


# Threshold computation
def mad(values: np.ndarray) -> float:
    """Median absolute deviation."""
    med = np.median(values)
    return float(np.median(np.abs(values - med)))


def compute_bin_thresholds(
    values: list[float],
    sided: str,
    n_mad: float,
) -> tuple[float, float, float | None, float] | None:
    """Compute median, MAD, and thresholds for one stratum.

    Parameters
    ----------
    values
        Raw metric values (NaN/inf are dropped).
    sided
        ``"one"`` for one-sided upper; ``"two"`` for two-sided.
    n_mad
        Threshold half-width in units of ``MAD_SCALE * MAD``.

    Returns
    -------
    tuple of (median, mad, lower, upper) or None
        ``lower`` is ``None`` for one-sided metrics.  Returns ``None`` if no
        finite values remain.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return None
    median = float(np.median(arr))
    m = mad(arr)
    half_width = n_mad * MAD_SCALE * m
    upper = median + half_width
    lower = (median - half_width) if sided == "two" else None
    return median, m, lower, upper


# Main
def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--butler-repo", required=True, help="Path or URI of the butler repository.")
    parser.add_argument("--collection", required=True, help="Butler collection to query.")
    parser.add_argument(
        "--pipeline",
        required=True,
        help="Path to the drp_pipe pipeline YAML, e.g. $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml.",
    )
    parser.add_argument(
        "--n-mad",
        type=float,
        default=DEFAULT_N_MAD,
        metavar="N",
        help="Threshold half-width in units of MAD_SCALE * MAD (default: %(default)s).",
    )
    parser.add_argument(
        "--tier1-fraction",
        type=float,
        default=DEFAULT_TIER1_FRACTION,
        help="Trip fraction above which a Tier 1 alarm is raised (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        help=(
            "Output CSV path.  Defaults to thresholds/<instrument>/<pipeline>.csv "
            "mirroring the pipeline path.  The PDF report is written alongside it."
        ),
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    # Metrics are bespoke to a given processing pipeline
    # so the output paths mirror the pipeline path.
    pipeline_path = Path(args.pipeline)
    instrument = pipeline_path.parent.name
    out_csv = (
        Path(args.output)
        if args.output
        else THRESHOLDS_DIR / instrument / pipeline_path.with_suffix(".csv").name
    )

    # Load supporting data
    metric_descs = load_metric_properties(METRIC_DESCRIPTIONS_PATH)
    bins = load_bins(instrument)

    gal_lat_edges: list[float] = bins["galactic_latitude"]["edges"]
    gal_lat_labels: list[str] = bins["galactic_latitude"]["labels"]

    # Parse pipeline using a temporary butler for graph resolution
    log.info("Parsing pipeline: %s", args.pipeline)
    test_dir = makeTestTempDir(str(_METRICINFO_DIR))
    try:
        with tempfile.TemporaryDirectory(dir=test_dir) as tmpdir:
            cfg = Config()
            cfg["registry", "db"] = f"sqlite:///{tmpdir}/gen3.sqlite3"
            tmp_repo = Butler.makeRepo(test_dir, cfg)
            tmp_butler = Butler.from_config(tmp_repo, writeable=True)
            table_info = walk_pipeline(args.pipeline, tmp_butler)
    finally:
        removeTestTempDir(test_dir)

    if not table_info:
        log.error("No metric tables found in pipeline graph. Nothing to do.")
        return 1

    # Check which tables actually exist in the collection
    log.info("Checking which metric tables exist in collection '%s'.", args.collection)
    butler = Butler(args.butler_repo, collections=args.collection)

    missing_tables: list[tuple[str, str]] = []
    existing_tables: list[str] = []
    for table_name, info in table_info.items():
        refs = list(butler.registry.queryDatasets(table_name, collections=args.collection))
        if refs:
            existing_tables.append(table_name)
        else:
            missing_tables.append((table_name, info["bundle_dataset_type"]))

    if missing_tables:
        log.warning(
            "%d metric table(s) are missing — run MakeMetricTableTask to generate them:\n%s",
            len(missing_tables),
            "\n".join(f"  table='{t}' (bundles='{b}')" for t, b in sorted(missing_tables)),
        )

    if not existing_tables:
        log.error(
            "No metric tables found in collection '%s'. "
            "Generate them with MakeMetricTableTask before running this script.",
            args.collection,
        )
        return 1

    log.info(
        "%d/%d metric table type(s) available; computing thresholds.",
        len(existing_tables),
        len(table_info),
    )


    # Load per-visit properties
    visit_properties = load_visit_properties(butler, args.collection)

    # Needed to look up metric descriptions (keyed by task_class) during
    # threshold computation, while the CSV uses the more readable task_name.
    task_name_to_class = {
        info["task_name"]: info["task_class"] for info in table_info.values()
    }

    # All values kept in memory: exact median/MAD needs the full distribution,
    # and the PDF report needs the raw values.  Revisit if memory becomes a bottleneck.
    data: dict[tuple, list[float]] = {}
    # Maps (task_name, atool, expanded_name) to template name, for description
    # lookups (metric_descs is keyed by template) and parquet output.
    metric_templates: dict[tuple[str, str, str], str] = {}
    skymap = None  # loaded lazily on first coadd-level table

    for table_name in existing_tables:
        info = table_info[table_name]
        task_name = info["task_name"]
        atools = info["atools"]

        refs = list(butler.registry.queryDatasets(table_name, collections=args.collection))
        if len(refs) > 1:
            log.warning(
                "Found %d datasets of type '%s'; expected 1.  This may indicate "
                "multiple processing runs in the same collection.  Using the last "
                "dataset — verify the collection contains only the intended run.",
                len(refs), table_name,
            )
        ref = refs[-1]
        log.info("Loading dataset of type '%s'.", table_name)

        try:
            table = butler.get(ref)
        except Exception as exc:
            log.warning("Could not load %s: %s", ref, exc)
            continue

        # Determine coadd vs. visit-level once per table
        is_coadd = "skymap" in ref.dataId.dimensions.names
        if is_coadd and skymap is None:
            skymap = butler.get("skyMap", skymap=ref.dataId["skymap"])

        # Pre-filter to metrics that have columns in this table.
        # Each entry is (atool, template, expanded_name, col): template is the
        # {band}-containing name from the pipeline config (used for description
        # lookups); expanded_name is the band-substituted form (used as the
        # metric identifier in the CSV and parquet).
        available: list[tuple[str, str, str, str]] = []
        for atool, metrics in atools.items():
            for template, bands in metrics.items():
                expansions = [template.format(band=b) for b in bands] if bands else [template]
                for expanded in expansions:
                    colname = f"{atool}_{expanded}"
                    if colname in table.colnames:
                        available.append((atool, template, expanded, colname))
                        metric_templates[(task_name, atool, expanded)] = template

        tract_gal_lat: dict[int, float] = {}  # store tract gal_lat centres within this table
        for row in table:
            if is_coadd:
                tract_id = int(row["tract"])
                if tract_id not in tract_gal_lat:
                    centre = skymap[tract_id].getCtrCoord()
                    tract_gal_lat[tract_id] = galactic_latitude(
                        centre.getRa().asDegrees(), centre.getDec().asDegrees()
                    )
                gal_lat = tract_gal_lat[tract_id]
            else:
                gal_lat = visit_properties.get(int(row["visit"]), {}).get("gal_lat")

            gal_lat_bin = (
                assign_bin(abs(gal_lat), gal_lat_edges, gal_lat_labels)
                if gal_lat is not None else "unknown"
            )
            # Traverse the columns of the row:
            for atool, template, expanded, col in available:
                key = (task_name, atool, expanded, gal_lat_bin)
                data.setdefault(key, []).append(float(row[col]))

    if not data:
        log.error("No metric values found — cannot compute thresholds.")
        return 1


    # Compute thresholds
    log.info(
        "Computing thresholds for %d (metric, stratum) combinations.", len(data)
    )
    rows: list[dict] = []
    now = datetime.now(timezone.utc).isoformat()

    for (
        task_name, atool, metric_name, gal_lat_bin
    ), values in sorted(data.items()):
        n = len(values)
        if n < MIN_VISITS_PER_STRATUM:
            log.warning(
                "%s / %s [b=%s]: only %d value(s) — thresholds may be unreliable.",
                atool, metric_name, gal_lat_bin, n,
            )

        task_class = task_name_to_class[task_name]
        template = metric_templates.get((task_name, atool, metric_name), metric_name)
        desc = metric_descs.get((task_class, atool, template), {})
        sided = desc.get("sided", "one")

        result = compute_bin_thresholds(values, sided, args.n_mad)
        if result is None:
            log.warning(
                "%s / %s [b=%s]: no finite values — skipping.",
                atool, metric_name, gal_lat_bin,
            )
            continue

        median, m, lower, upper = result

        if m == 0.0:
            log.warning(
                "%s / %s [b=%s]: MAD=0 — all values identical?",
                atool, metric_name, gal_lat_bin,
            )

        rows.append({
            "task_name": task_name,
            "atool": atool,
            "metric_name": metric_name,
            "gal_lat_bin": gal_lat_bin,
            "lower": lower,
            "upper": upper,
            "median": median,
            "mad": m,
            "n_values": n,
            "tier1_fraction": args.tier1_fraction,
            "computed_at": now,
        })


    # Validate output
    df = pd.DataFrame(rows)

    n_nan_upper = int(df["upper"].isna().sum())
    if n_nan_upper:
        log.error("%d rows have a NaN upper threshold.", n_nan_upper)

    two_sided_mask = df["lower"].notna()
    if two_sided_mask.any():
        n_nan_lower = int(df.loc[two_sided_mask, "lower"].isna().sum())
        if n_nan_lower:
            log.error("%d two-sided rows have a NaN lower threshold.", n_nan_lower)

    n_zero_mad = int((df["mad"] == 0.0).sum())
    if n_zero_mad:
        log.warning("%d rows have MAD=0.", n_zero_mad)

    n_unknown = int((df["gal_lat_bin"] == "unknown").sum())
    if n_unknown:
        log.warning(
            "%d rows have an unknown stratum — conditions data may be incomplete.",
            n_unknown,
        )

    log.info(
        "Validation complete: %d rows total, %d two-sided, %d one-sided.",
        len(df),
        int(two_sided_mask.sum()),
        int((~two_sided_mask).sum()),
    )


    # Write CSV, parquet + PDF.
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    log.info("Thresholds written to %s.", out_csv)

    # Save raw values as long-format parquet in the current working directory
    pipeline_stem = pipeline_path.stem
    out_parquet = Path.cwd() / f"{pipeline_stem}_values.parquet"
    values_records = []
    for (task_name, atool, metric_name, gal_lat_bin), vals in data.items():
        task_class = task_name_to_class[task_name]
        template = metric_templates.get((task_name, atool, metric_name), metric_name)
        for v in vals:
            values_records.append({
                "task_name": task_name,
                "task_class": task_class,
                "atool": atool,
                "metric_name": metric_name,
                "template_metric_name": template,  # {band} form; used by generate_report for description lookups
                "gal_lat_bin": gal_lat_bin,
                "value": v,
            })
    pd.DataFrame(values_records).to_parquet(out_parquet, index=False)
    log.info("Raw values written to %s.", out_parquet)

    # Generate PDF report in the current working directory
    out_pdf = Path.cwd() / f"{pipeline_stem}.pdf"
    generate_report(
        thresholds_path=out_csv,
        values_path=out_parquet,
        output_path=out_pdf,
        pipeline_name=f"{instrument}/{pipeline_path.name}",
        collection=args.collection,
    )

    if missing_tables:
        log.warning(
            "%d metric table type(s) were missing and excluded. "
            "Re-run after generating them to include all metrics.",
            len(missing_tables),
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
