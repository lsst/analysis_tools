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
"""Generate a threshold report PDF from a threshold CSV and values parquet.

Can be run standalone to regenerate the PDF without recomputing thresholds,
or called from compute_thresholds.py after threshold computation.

The values parquet is a long-format table with columns:
  task_class, atool, metric_name, gal_lat_bin, value

Examples
--------
::

    python generate_report.py \\
        --thresholds thresholds/LSSTCam/DRP.csv \\
        --values DRP_values.parquet \\
        --output DRP.pdf
"""
from __future__ import annotations

import argparse
import logging
import math
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd
import yaml

log = logging.getLogger(__name__)

_SCRIPTS_DIR = Path(__file__).resolve().parent
METRIC_DESCRIPTIONS_PATH = _SCRIPTS_DIR.parent / "metricDescriptions.yaml"

METRICS_PER_PAGE = 4
DESCRIPTION_WRAP_WIDTH = 38
LEVEL_TYPES = {"visit", "coadd"}
BAND_ORDER = {"u": 0, "g": 1, "r": 2, "i": 3, "z": 4, "y": 5}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def load_metric_descriptions(path: Path) -> dict[tuple[str, str, str], dict]:
    """Load and flatten metricDescriptions.yaml.

    Parameters
    ----------
    path
        Path to metricDescriptions.yaml.

    Returns
    -------
    dict
        Keyed by ``(task_class, atool, metric_name)``.
    """
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    result: dict[tuple[str, str, str], dict] = {}
    for task_class, task_body in raw["tasks"].items():
        for atool, atool_body in task_body["atools"].items():
            for metric_name, props in atool_body.items():
                result[(task_class, atool, metric_name)] = props
    return result


def _metric_type_label(metric_types) -> str:
    """Return a category label for a metric, excluding level types."""
    if metric_types is None:
        return "other"
    if isinstance(metric_types, str):
        types = [metric_types]
    else:
        types = list(metric_types)
    categories = [t for t in types if t not in LEVEL_TYPES]
    return ", ".join(sorted(categories)) if categories else "other"


def _metric_level(metric_types) -> str:
    """Return 'visit', 'coadd', or 'misc' for a metric based on its metricTypes."""
    if metric_types is None:
        return "misc"
    types = [metric_types] if isinstance(metric_types, str) else list(metric_types)
    if "coadd" in types:
        return "coadd"
    if "visit" in types:
        return "visit"
    return "misc"


def _stratum_label(gal: str) -> str:
    return f"b={gal}"


def _metric_sort_key(key: tuple[str, str, str]) -> tuple:
    """Sort key that groups band-prefixed metrics by base name, ordered u,g,r,i,z."""
    task_name, atool, metric_name = key
    parts = metric_name.split("_", 1)
    if len(parts) == 2 and parts[0] in BAND_ORDER:
        band, base = parts
        return (task_name, atool, base, BAND_ORDER[band], metric_name)
    return (task_name, atool, metric_name, -1, metric_name)


# ---------------------------------------------------------------------------
# Page rendering
# ---------------------------------------------------------------------------


def _render_page(
    pdf: PdfPages,
    page_metrics: list[tuple[str, str, str]],
    strata: list[tuple[str,]],
    values_df: pd.DataFrame,
    thresh_df: pd.DataFrame,
    name_descs: dict[tuple[str, str, str], dict],
    level: str,
    type_label: str,
    page_label: str,
    pipeline_name: str,
) -> None:
    n_rows = len(page_metrics)
    n_strata = len(strata)
    colours = plt.cm.tab10.colors[:n_strata]

    col_widths = [2.4, 1.8 * n_strata]
    fig_width = sum(col_widths) + 0.4
    fig_height = 2.5 * n_rows + 0.5

    fig = plt.figure(figsize=(fig_width, fig_height))
    level_label = "unknown" if level == "misc" else level
    fig.suptitle(
        f"Level: {level_label}  |  Type: {type_label}{page_label}  |  {pipeline_name}",
        fontsize=10, fontweight="bold", y=0.98,
    )

    outer_gs = gridspec.GridSpec(
        n_rows, 2,
        figure=fig,
        width_ratios=col_widths,
        top=0.93,
        hspace=0.65,
        wspace=0.08,
    )

    row_info: list[tuple] = []

    for row_idx, (task_name, atool, metric_name) in enumerate(page_metrics):
        desc_props = name_descs.get((task_name, atool, metric_name), {})
        units = desc_props.get("units", "")
        description = desc_props.get("description", "No description available.")
        metric_types = desc_props.get("metricTypes", [])
        if isinstance(metric_types, str):
            metric_types = [metric_types]

        # ---- Info panel ----
        info_ax = fig.add_subplot(outer_gs[row_idx, 0])
        inner_gs = gridspec.GridSpecFromSubplotSpec(
            1, n_strata, subplot_spec=outer_gs[row_idx, 1], wspace=0,
        )
        info_ax.axis("off")

        wrapped_desc = "\n".join(textwrap.wrap(description, width=DESCRIPTION_WRAP_WIDTH))
        meta_lines = []
        if units:
            meta_lines.append(("Units:", units))
        meta_lines.append(("Sided:", desc_props.get("sided", "one")))
        if metric_types:
            meta_lines.append(("Type:", ", ".join(metric_types)))
        meta_lines.append(("Task:", task_name))
        meta_lines.append(("ATool:", atool))

        info_ax.text(
            0.0, 1.0, metric_name,
            transform=info_ax.transAxes,
            fontsize=7, fontweight="bold",
            va="top", ha="left", clip_on=False,
        )
        info_ax.text(
            0.0, 0.82, wrapped_desc,
            transform=info_ax.transAxes,
            fontsize=6, va="top", ha="left",
            clip_on=False, linespacing=1.4,
        )

        meta_y = 0.18
        line_step = 0.13
        for label, value in meta_lines:
            info_ax.text(
                0.0, meta_y, label,
                transform=info_ax.transAxes,
                fontsize=6, fontweight="bold", va="top", ha="left", clip_on=False,
            )
            info_ax.text(
                0.28, meta_y, value,
                transform=info_ax.transAxes,
                fontsize=6, va="top", ha="left", clip_on=False,
            )
            meta_y -= line_step

        row_info.append((info_ax, meta_y - 0.05))

        # ---- Histogram panels ----
        hist_axes = []
        for col_idx, (gal_bin,) in enumerate(strata):
            share_ax = hist_axes[0] if col_idx > 0 else None
            ax = fig.add_subplot(inner_gs[0, col_idx], sharey=share_ax, sharex=share_ax)
            hist_axes.append(ax)

            colour = colours[col_idx]

            # Raw values for this metric + stratum
            mask = (
                (values_df["task_name"] == task_name)
                & (values_df["atool"] == atool)
                & (values_df["metric_name"] == metric_name)
                & (values_df["gal_lat_bin"] == gal_bin)
            )
            vals = values_df.loc[mask, "value"].dropna().values

            if len(vals) > 0:
                ax.hist(
                    vals,
                    bins=min(30, max(5, len(vals) // 5)),
                    color=colour,
                    alpha=0.6,
                    edgecolor="none",
                    density=True,
                )

            # Thresholds for this metric + stratum
            tmask = (
                (thresh_df["task_name"] == task_name)
                & (thresh_df["atool"] == atool)
                & (thresh_df["metric_name"] == metric_name)
                & (thresh_df["gal_lat_bin"] == gal_bin)
            )
            trows = thresh_df.loc[tmask]
            if not trows.empty:
                row = trows.iloc[0]
                ax.axvline(row["median"], color=colour, lw=1.2, ls="-")
                ax.axvline(row["upper"], color=colour, lw=1.0, ls="--")
                if pd.notna(row["lower"]):
                    ax.axvline(row["lower"], color=colour, lw=1.0, ls="--")

            ax.tick_params(labelsize=5)
            if units:
                ax.set_xlabel(units, fontsize=5, labelpad=2)
            ax.text(
                0.98, 0.97, f"n={len(vals)}",
                transform=ax.transAxes,
                fontsize=5, ha="right", va="top",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5),
            )

            if row_idx == 0:
                ax.set_title(_stratum_label(gal_bin,), fontsize=6)

            if col_idx > 0:
                ax.tick_params(left=False, labelleft=False)

    # Divider lines between rows
    for info_ax, divider_y in row_info[:-1]:
        trans = blended_transform_factory(fig.transFigure, info_ax.transAxes)
        line = Line2D(
            [0.1, 0.9], [divider_y, divider_y],
            transform=trans,
            color="lightgrey",
            linewidth=0.8,
            clip_on=False,
        )
        fig.add_artist(line)

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_report(
    thresholds_path: Path,
    values_path: Path,
    output_path: Path,
    pipeline_name: str = "",
    collection: str = "",
    test: bool = False,
) -> None:
    """Generate a threshold report PDF.

    Parameters
    ----------
    thresholds_path
        Path to the threshold CSV produced by compute_thresholds.py.
    values_path
        Path to the long-format parquet of raw metric values.
    output_path
        Output PDF path.
    pipeline_name
        Pipeline name shown in page headers (optional).
    collection
        Butler collection shown on the title page (optional).
    test
        If True, render only the first five pages of the report for testing purposes.
    """
    metric_descs = load_metric_descriptions(METRIC_DESCRIPTIONS_PATH)
    thresh_df = pd.read_csv(thresholds_path)
    values_df = pd.read_parquet(values_path)

    # Build a description lookup keyed by (task_name, atool, expanded_metric_name).
    # metricDescriptions.yaml is keyed by (task_class, atool, template_metric_name),
    # where template_metric_name may contain {band} (e.g. "{band}_AM1").  The parquet
    # stores both the expanded name (e.g. "r_AM1") and the template, so we can look
    # up the description directly without needing to infer the band substitution.
    metric_lookup = (
        values_df[["task_name", "task_class", "atool", "metric_name", "template_metric_name"]]
        .drop_duplicates()
        .itertuples(index=False)
    )
    name_descs = {
        (task_name, atool, metric_name): metric_descs.get((task_class, atool, template), {})
        for task_name, task_class, atool, metric_name, template in metric_lookup
    }

    # Determine the ordered list of strata from the threshold CSV
    strata: list[tuple[str,]] = (
        thresh_df[["gal_lat_bin",]]
        .drop_duplicates()
        .apply(tuple, axis=1)
        .tolist()
    )
    strata = sorted(set(strata))

    # Group metrics by (level, category type); visit first, then coadd, then misc.
    _LEVEL_ORDER = {"visit": 0, "coadd": 1, "misc": 2}
    from collections import defaultdict
    type_to_metrics: dict[tuple[str, str], list[tuple]] = defaultdict(list)
    seen: set[tuple] = set()
    for _, row in thresh_df.iterrows():
        key = (row["task_name"], row["atool"], row["metric_name"])
        if key in seen:
            continue
        seen.add(key)
        desc = name_descs.get(key, {})
        metric_types = desc.get("metricTypes")
        label = _metric_type_label(metric_types)
        level = _metric_level(metric_types)
        type_to_metrics[(level, label)].append(key)

    with PdfPages(output_path) as pdf:
        # Title page
        fig, ax = plt.subplots(figsize=(8.27, 11.69))
        ax.axis("off")
        ax.text(0.5, 0.65, "Metric Threshold Report",
                ha="center", va="center", fontsize=20, fontweight="bold",
                transform=ax.transAxes)
        ax.text(0.5, 0.55, f"Pipeline: {pipeline_name}",
                ha="center", va="center", fontsize=12, transform=ax.transAxes)
        ax.text(0.5, 0.50, f"Collection: {collection}",
                ha="center", va="center", fontsize=10, transform=ax.transAxes)
        ax.text(0.5, 0.44,
                datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
                ha="center", va="center", fontsize=9, color="grey",
                transform=ax.transAxes)
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        # Metric pages grouped by metricType, visit-level first then coadd then misc
        pages_rendered = 0
        for level, type_label in sorted(type_to_metrics, key=lambda x: (_LEVEL_ORDER[x[0]], x[1])):
            if test and pages_rendered >= 5:
                break
            metrics = sorted(type_to_metrics[(level, type_label)], key=_metric_sort_key)
            n_pages = math.ceil(len(metrics) / METRICS_PER_PAGE)
            for page_idx in range(n_pages):
                if test and pages_rendered >= 5:
                    break
                page_metrics = metrics[
                    page_idx * METRICS_PER_PAGE: (page_idx + 1) * METRICS_PER_PAGE
                ]
                page_label = f"  (page {page_idx + 1}/{n_pages})" if n_pages > 1 else ""
                _render_page(
                    pdf, page_metrics, strata,
                    values_df, thresh_df, name_descs,
                    level, type_label, page_label, pipeline_name,
                )
                pages_rendered += 1

    log.info("PDF report written to %s.", output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--thresholds", required=True,
        help="Path to threshold CSV (e.g. metricInfo/thresholds/LSSTCam/DRP.csv).",
    )
    parser.add_argument(
        "--values", required=True,
        help="Path to values parquet (e.g. ./DRP_values.parquet).",
    )
    parser.add_argument(
        "--output",
        help="Output PDF path (default: <values stem>.pdf in current directory).",
    )
    parser.add_argument("--pipeline", default="", help="Pipeline name for report header.")
    parser.add_argument("--collection", default="", help="Butler collection for report header.")
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    parser.add_argument(
        "--dummy-report", action="store_true",
        help="Render only the first five pages for testing purposes.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    values_path = Path(args.values)
    out_path = (
        Path(args.output) if args.output
        else Path.cwd() / values_path.with_suffix(".pdf").name
    )

    generate_report(
        thresholds_path=Path(args.thresholds),
        values_path=values_path,
        output_path=out_path,
        pipeline_name=args.pipeline,
        collection=args.collection,
        test=args.dummy_report,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
