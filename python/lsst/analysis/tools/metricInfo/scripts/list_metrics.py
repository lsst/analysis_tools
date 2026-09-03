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
"""List all metrics produced by a drp_pipe pipeline, with their descriptions.

Parses a drp_pipe pipeline YAML to discover metric-producing tasks (the same
graph walk used by ``compute_thresholds.py``), then joins the resulting
``(task_class, atool, metric_name)`` tuples against ``metricDescriptions.yaml``
to report each metric's units, sidedness, and metricTypes — and against
``primaryMetrics.yaml`` to flag which ones are in the curated Tier 1 set.

Pipeline graph resolution needs a butler registry, so a throwaway sqlite one
is created in a temp directory purely to satisfy that — no real butler repo
or collection is contacted, and no data is read. Metric name templates are
left as ``{band}`` placeholders; each template is listed once, alongside the
bands the task expands it to.

Examples
--------
::

    python list_metrics.py --pipeline $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml
    python list_metrics.py --pipeline $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml \\
        --output DRP_metrics.csv
    python list_metrics.py --pipeline $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml \\
        --html DRP_metrics.html
"""
from __future__ import annotations

import argparse
import json
import logging
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

from compute_thresholds import (
    METRIC_DESCRIPTIONS_PATH,
    extract_task_metrics,
    load_metric_properties,
)

from lsst.daf.butler import Butler, Config
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.pipe.base import Pipeline

log = logging.getLogger(__name__)

_METRICINFO_DIR = Path(__file__).resolve().parent.parent
PRIMARY_METRICS_PATH = _METRICINFO_DIR / "primaryMetrics.yaml"


def load_primary_metrics(path: Path) -> set[tuple[str, str, str]]:
    """Load primaryMetrics.yaml into a set of identifying tuples.

    Unlike metricDescriptions.yaml, this file carries no per-metric
    properties — just membership — so each atool maps to a plain list of
    metric names rather than a name -> properties dict.

    Parameters
    ----------
    path
        Path to primaryMetrics.yaml.

    Returns
    -------
    set
        ``(task_class, atool, metric_name)`` tuples in the primary set.
    """
    if not path.exists():
        return set()
    with open(path) as file:
        md = yaml.safe_load(file)
    return {
        (task_class, atool, metric_name)
        for task_class, task_cfg in md["tasks"].items()
        for atool, metric_names in task_cfg["atools"].items()
        for metric_name in metric_names
    }


def walk_pipeline_all(pipeline_path: str, butler: Butler) -> dict[str, dict]:
    """Walk a pipeline graph and return every metric-producing task.

    Unlike `compute_thresholds.walk_pipeline`, this doesn't require the
    metric bundle to reach a `MakeMetricTableTask` — that requirement exists
    there because threshold computation needs real, tabulated data, which is
    irrelevant here: this never touches real data, only the pipeline graph
    and metricDescriptions.yaml. A metric bundle with no table still has
    every property metricDescriptions.yaml can give it; it just can't
    currently feed threshold computation, which each entry records via
    ``has_table``.

    Parameters
    ----------
    pipeline_path
        Path to the drp_pipe pipeline YAML.
    butler
        Butler instance used for pipeline graph resolution.

    Returns
    -------
    dict
        Keyed by task_name (unique per pipeline). Each value is a dict with
        keys ``task_name``, ``task_class``, ``atools``, ``has_table``.
    """
    pipeline = Pipeline.from_uri(pipeline_path)
    pipeline_graph = pipeline.to_graph(registry=butler.registry, visualization_only=True)

    bundle_to_table: dict[str, str] = {}
    for task in pipeline_graph.tasks.values():
        if "MakeMetricTableTask" not in task.task_class_name:
            continue
        bundle_to_table[task.inputs["data"].dataset_type_name] = task.outputs["metricTable"].dataset_type_name

    result: dict[str, dict] = {}
    for task_name, task in pipeline_graph.tasks.items():
        if "metrics" not in task.outputs:
            continue
        bundle_name = task.outputs["metrics"].dataset_type_name
        atools = extract_task_metrics(task.config.toDict())
        if not atools:
            continue
        result[task_name] = {
            "task_name": task_name,
            "task_class": task.task_class_name,
            "atools": atools,
            "has_table": bundle_name in bundle_to_table,
        }
    return result


def build_metric_table(pipeline_path: str) -> pd.DataFrame:
    """Build a table of all metrics produced by a pipeline, with descriptions.

    Parameters
    ----------
    pipeline_path
        Path to the drp_pipe pipeline YAML.

    Returns
    -------
    pandas.DataFrame
        One row per ``(task_name, atool, metric_name)``, with columns
        ``task_name``, ``task_class``, ``atool``, ``metric_name``, ``bands``,
        ``units``, ``sided``, ``metricTypes``, ``description``, ``is_primary``,
        ``has_table``. Metric names containing ``{band}`` are listed once as
        a template, not expanded.
    """
    metric_descs = load_metric_properties(METRIC_DESCRIPTIONS_PATH)
    primary_metrics = load_primary_metrics(PRIMARY_METRICS_PATH)

    # Graph resolution requires a registry, but not real data: a throwaway
    # sqlite one, discarded immediately, stands in for it.
    test_dir = makeTestTempDir(str(_METRICINFO_DIR))
    try:
        with tempfile.TemporaryDirectory(dir=test_dir) as tmpdir:
            cfg = Config()
            cfg["registry", "db"] = f"sqlite:///{tmpdir}/gen3.sqlite3"
            tmp_repo = Butler.makeRepo(test_dir, cfg)
            tmp_butler = Butler.from_config(tmp_repo, writeable=True)
            task_info = walk_pipeline_all(pipeline_path, tmp_butler)
    finally:
        removeTestTempDir(test_dir)

    rows = []
    for info in task_info.values():
        task_name = info["task_name"]
        task_class = info["task_class"]
        has_table = info["has_table"]
        for atool, metrics in info["atools"].items():
            for metric_name, bands in metrics.items():
                props = metric_descs.get((task_class, atool, metric_name), {})
                if not props:
                    log.warning(
                        "No metricDescriptions.yaml entry for (%s, %s, %s) — "
                        "leaving description/units/sided/metricTypes blank.",
                        task_class, atool, metric_name,
                    )
                # `.get(key, default)` only substitutes when the key is absent — several
                # metricDescriptions.yaml entries spell "no units" as an explicit `units:
                # null` rather than omitting the key, which `.get` happily hands back as
                # None. `or` catches that too, alongside the genuinely-missing case.
                metric_types = props.get("metricTypes") or []
                if isinstance(metric_types, str):
                    metric_types = [metric_types]
                rows.append(
                    {
                        "task_name": task_name,
                        "task_class": task_class,
                        "atool": atool,
                        "metric_name": metric_name,
                        "bands": ",".join(bands),
                        "units": props.get("units") or "",
                        "sided": props.get("sided") or "",
                        "metricTypes": ",".join(metric_types),
                        "description": props.get("description", ""),
                        "is_primary": (task_class, atool, metric_name) in primary_metrics,
                        "has_table": has_table,
                    }
                )

    df = pd.DataFrame(rows).sort_values(["task_name", "atool", "metric_name"])
    return df.reset_index(drop=True)


# Approximate ugrizy filter-band colours, used only to give the "bands" column
# a quick visual identity — not a claim of any official Rubin palette.
_BAND_COLORS = {
    "u": "#3f8fd1",
    "g": "#2f9e63",
    "r": "#d1552f",
    "i": "#9c3a3a",
    "z": "#7a4fb5",
    "y": "#232323",
}

_HTML_TEMPLATE = """<!doctype html>
<title>__TITLE__</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,600;8..60,700&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap');

:root {
  --bg: #f5f6f9;
  --surface: #ffffff;
  --surface-2: #eceef4;
  --ink: #171b26;
  --ink-dim: #5c6478;
  --line: #dde1ea;
  --accent: #b5591f;
  --accent-ink: #ffffff;
  --good: #2f7d5c;
  --gap: #a83e3e;
  --focus: #3f6fd1;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #10131b;
    --surface: #171b26;
    --surface-2: #1e2330;
    --ink: #e7e9ef;
    --ink-dim: #98a0b5;
    --line: #2a3040;
    --accent: #e2903f;
    --accent-ink: #171307;
    --good: #4fae83;
    --gap: #e2726a;
    --focus: #6f96ea;
  }
}
:root[data-theme="dark"] {
  --bg: #10131b;
  --surface: #171b26;
  --surface-2: #1e2330;
  --ink: #e7e9ef;
  --ink-dim: #98a0b5;
  --line: #2a3040;
  --accent: #e2903f;
  --accent-ink: #171307;
  --good: #4fae83;
  --gap: #e2726a;
  --focus: #6f96ea;
}

* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--ink);
  font: 15px/1.5 "IBM Plex Sans", system-ui, sans-serif;
}
::selection { background: var(--accent); color: var(--accent-ink); }
:focus-visible { outline: 2px solid var(--focus); outline-offset: 2px; }

header {
  max-width: 78rem;
  margin: 0 auto;
  padding: 2.75rem 1.75rem 1.25rem;
}
.eyebrow {
  font: 500 12px/1 "IBM Plex Mono", monospace;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--accent);
}
h1 {
  font: 700 2.35rem/1.1 "Source Serif 4", Georgia, serif;
  text-wrap: balance;
  margin: 0.4rem 0 0.6rem;
}
.stats {
  color: var(--ink-dim);
  font-variant-numeric: tabular-nums;
  font-size: 0.95rem;
}
.stats strong { color: var(--ink); font-weight: 600; }
.stats .good { color: var(--good); font-weight: 600; }
.stats .gap { color: var(--gap); font-weight: 600; }
.stats .primary-count { color: var(--accent); font-weight: 600; }

.toolbar {
  position: sticky;
  top: 0;
  z-index: 5;
  background: var(--bg);
  border-bottom: 1px solid var(--line);
  padding: 0.85rem 1.75rem;
}
.toolbar-inner {
  max-width: 78rem;
  margin: 0 auto;
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 0.6rem;
}
#search {
  flex: 1 1 16rem;
  min-width: 12rem;
  padding: 0.5rem 0.8rem;
  border: 1px solid var(--line);
  border-radius: 7px;
  background: var(--surface);
  color: var(--ink);
  font: 500 0.9rem "IBM Plex Mono", monospace;
}
#search::placeholder { color: var(--ink-dim); }

.chip {
  border: 1px solid var(--line);
  background: var(--surface);
  color: var(--ink-dim);
  border-radius: 999px;
  padding: 0.32rem 0.75rem;
  font: 500 0.78rem "IBM Plex Sans", sans-serif;
  letter-spacing: 0.01em;
  cursor: pointer;
  user-select: none;
  white-space: nowrap;
}
.chip[aria-pressed="true"] {
  background: var(--accent);
  border-color: var(--accent);
  color: var(--accent-ink);
}
.chip .n { opacity: 0.7; font-variant-numeric: tabular-nums; }
.chip-group { display: flex; flex-wrap: wrap; gap: 0.4rem; }
.divider { width: 1px; align-self: stretch; background: var(--line); margin: 0 0.2rem; }

main { max-width: 78rem; margin: 0 auto; padding: 1.25rem 1.75rem 4rem; }

/* `overflow-x: auto` alone would implicitly compute overflow-y to `auto`
   too (a lone `visible` axis isn't allowed once the other isn't visible),
   silently turning this into its own scroll container — which then traps
   position:sticky, since sticky sticks to the nearest scrolling ancestor,
   and this box never scrolls the header past that ancestor's own bounds.
   Making that bounded scroll deliberate — a fixed max-height, scrolling
   in both axes — gives the header a real scrollport to stick within. */
.table-wrap {
  overflow: auto;
  max-height: calc(100vh - 16rem);
  border: 1px solid var(--line);
  border-radius: 10px;
  background: var(--surface);
}
table { width: 100%; border-collapse: collapse; min-width: 62rem; }
thead th {
  position: sticky;
  top: 0;
  background: var(--surface-2);
  text-align: left;
  font: 600 0.72rem "IBM Plex Sans", sans-serif;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  color: var(--ink-dim);
  padding: 0.65rem 0.9rem;
  border-bottom: 1px solid var(--line);
  cursor: pointer;
  white-space: nowrap;
}
thead th:hover { color: var(--ink); }
thead th .arrow { opacity: 0.4; margin-left: 0.15rem; }
thead th[aria-sort="ascending"] .arrow,
thead th[aria-sort="descending"] .arrow { opacity: 1; color: var(--accent); }

tbody td {
  padding: 0.65rem 0.9rem;
  border-bottom: 1px solid var(--line);
  vertical-align: top;
  font-size: 0.86rem;
}
tbody tr:last-child td { border-bottom: none; }
tbody tr:hover { background: var(--surface-2); }

.col-cov { width: 1.6rem; text-align: center; }
.dot { display: inline-block; width: 8px; height: 8px; border-radius: 999px; margin-top: 0.35rem; }
.dot.good { background: var(--good); }
.dot.gap { background: var(--gap); }
.star { color: var(--accent); font-size: 0.9rem; }
.table-icon { display: block; margin: 0.15rem auto 0; }
.table-icon.good { color: var(--good); }
.table-icon.gap { color: var(--gap); }

.metric-name { font: 500 0.86rem "IBM Plex Mono", monospace; color: var(--ink); }

.task-name { font-size: 0.82rem; }
.task-class { color: var(--ink-dim); font-size: 0.72rem; margin-top: 0.15rem; word-break: break-word; }
.task-atool { color: var(--ink-dim); font-size: 0.78rem; margin-top: 0.15rem; }

.bands { display: flex; gap: 0.3rem; flex-wrap: wrap; }
.band-dot { width: 9px; height: 9px; border-radius: 999px; border: 1px solid var(--line); }

.tags { display: flex; gap: 0.3rem; flex-wrap: wrap; }
.tag {
  font: 500 0.7rem "IBM Plex Mono", monospace;
  padding: 0.1rem 0.4rem;
  border-radius: 5px;
  background: var(--surface-2);
  color: var(--ink-dim);
  border: 1px solid var(--line);
}
.tag.level { color: var(--accent); border-color: var(--accent); }

.units-sided { font: 400 0.78rem "IBM Plex Mono", monospace; color: var(--ink-dim); white-space: nowrap; }

.description { max-width: 30ch; color: var(--ink); }
.description.empty { color: var(--gap); font-style: italic; }
.description.empty::before { content: "no catalog entry yet"; }

.count-line { color: var(--ink-dim); font-size: 0.82rem; margin: 0.75rem 0.2rem; }
.empty-state { padding: 3rem 1rem; text-align: center; color: var(--ink-dim); }
[hidden] { display: none !important; }
</style>

<header>
  <div class="eyebrow">__PIPELINE_LABEL__ · analysis_tools metricInfo</div>
  <h1>__H1__</h1>
  <p class="stats">
    <strong>__TOTAL__</strong> metrics discovered ·
    <span class="good">__DOCUMENTED__ documented</span> ·
    <span class="gap">__GAPS__ gaps</span>
    (__COVERAGE_PCT__% coverage) ·
    <span class="primary-count">__PRIMARY_COUNT__ primary</span> ·
    <span class="gap">__NO_TABLE__ not tabulated</span> · generated __GENERATED__
  </p>
</header>

<div class="toolbar">
  <div class="toolbar-inner">
    <input id="search" type="search" placeholder="Search metric, task, atool, description…" autocomplete="off">
    <div class="divider"></div>
    <div class="chip-group" id="level-chips"></div>
    <div class="divider"></div>
    <div class="chip-group" id="type-chips"></div>
    <div class="divider"></div>
    <div class="chip-group">
      <button class="chip" id="gaps-only" aria-pressed="false" type="button">Gaps only</button>
      <button class="chip" id="primary-only" aria-pressed="false" type="button">Primary only</button>
      <button class="chip" id="no-table-only" aria-pressed="false" type="button">Not tabulated</button>
    </div>
  </div>
</div>

<main>
  <p class="count-line" id="count-line"></p>
  <div class="table-wrap">
    <table>
      <thead>
        <tr>
          <th class="col-cov" data-key="_covered" title="Has a metricDescriptions.yaml entry"></th>
          <th class="col-cov" data-key="is_primary" title="In the curated Tier 1 primary-metric set"></th>
          <th class="col-cov" data-key="has_table" title="Fed into a MakeMetricTableTask — can currently have thresholds computed"></th>
          <th data-key="metric_name">Metric <span class="arrow">↕</span></th>
          <th data-key="task_name">Pipetask / Task / ATool <span class="arrow">↕</span></th>
          <th data-key="bands">Bands <span class="arrow">↕</span></th>
          <th data-key="metricTypes">Type <span class="arrow">↕</span></th>
          <th data-key="units">Units / sided <span class="arrow">↕</span></th>
          <th data-key="description">Description <span class="arrow">↕</span></th>
        </tr>
      </thead>
      <tbody id="rows"></tbody>
    </table>
  </div>
  <div class="empty-state" id="empty-state" hidden>No metrics match the current filters.</div>
</main>

<script>
const DATA = __DATA_JSON__;
const BAND_COLORS = __BAND_COLORS_JSON__;

// Every task class lives under this same package; showing it on each row
// is just noise.
const TASK_CLASS_PREFIX = "lsst.analysis.tools.tasks.";
for (const row of DATA) {
  if (row.task_class.startsWith(TASK_CLASS_PREFIX)) row.task_class = row.task_class.slice(TASK_CLASS_PREFIX.length);
}

const state = {
  search: "", level: null, types: new Set(), gapsOnly: false, primaryOnly: false, noTableOnly: false,
  sortKey: "task_name", sortDir: 1,
};

function levelOf(row) {
  const types = row.metricTypes.split(",").filter(Boolean);
  if (types.includes("visit")) return "visit";
  if (types.includes("coadd")) return "coadd";
  return null;
}
function categoryTypes(row) {
  return row.metricTypes.split(",").filter(t => t && t !== "visit" && t !== "coadd");
}

const levels = ["visit", "coadd"].filter(lv => DATA.some(r => levelOf(r) === lv));
const allTypes = [...new Set(DATA.flatMap(categoryTypes))].sort();

function buildChips(container, values, isActive, onToggle, countFn) {
  container.innerHTML = "";
  for (const v of values) {
    const btn = document.createElement("button");
    btn.className = "chip";
    btn.type = "button";
    btn.setAttribute("aria-pressed", isActive(v));
    btn.innerHTML = v + ' <span class="n">' + countFn(v) + "</span>";
    btn.addEventListener("click", () => { onToggle(v); render(); });
    container.appendChild(btn);
  }
}

function refreshChips() {
  buildChips(
    document.getElementById("level-chips"), levels,
    v => state.level === v,
    v => { state.level = state.level === v ? null : v; },
    v => DATA.filter(r => levelOf(r) === v).length
  );
  buildChips(
    document.getElementById("type-chips"), allTypes,
    v => state.types.has(v),
    v => { state.types.has(v) ? state.types.delete(v) : state.types.add(v); },
    v => DATA.filter(r => categoryTypes(r).includes(v)).length
  );
}

function bandDots(bandsStr) {
  const bands = bandsStr.split(",").filter(Boolean);
  if (!bands.length) return "";
  return '<span class="bands">' + bands.map(b =>
    '<span class="band-dot" title="' + b + '" style="background:' + (BAND_COLORS[b] || "#888") + '"></span>'
  ).join("") + "</span>";
}

function typeTags(metricTypes) {
  const types = metricTypes.split(",").filter(Boolean);
  return '<span class="tags">' + types.map(t =>
    '<span class="tag' + ((t === "visit" || t === "coadd") ? " level" : "") + '">' + t + "</span>"
  ).join("") + "</span>";
}

function escapeHtml(s) {
  return String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

function tableIcon(hasTable) {
  const cls = hasTable ? "good" : "gap";
  const title = hasTable
    ? "Fed into a MakeMetricTableTask — can currently have thresholds computed"
    : "Not tabulated — no MakeMetricTableTask consumes this metric yet, so it can't currently have thresholds computed";
  return '<svg class="table-icon ' + cls + '" viewBox="0 0 12 12" width="11" height="11" aria-hidden="true">' +
    '<title>' + title + '</title>' +
    '<rect x="1" y="1" width="10" height="10" rx="1.5" fill="none" stroke="currentColor" stroke-width="1.2"/>' +
    '<line x1="1" y1="6" x2="11" y2="6" stroke="currentColor" stroke-width="1"/>' +
    '<line x1="6" y1="1" x2="6" y2="11" stroke="currentColor" stroke-width="1"/>' +
    '</svg>';
}

function matches(row) {
  if (state.level && levelOf(row) !== state.level) return false;
  if (state.types.size && ![...state.types].every(t => categoryTypes(row).includes(t))) return false;
  if (state.gapsOnly && row.description) return false;
  if (state.primaryOnly && !row.is_primary) return false;
  if (state.noTableOnly && row.has_table) return false;
  if (state.search) {
    const hay = (row.metric_name + " " + row.task_name + " " + row.atool + " " + row.description).toLowerCase();
    if (!hay.includes(state.search)) return false;
  }
  return true;
}

function sortValue(row) {
  const v = row[state.sortKey];
  return typeof v === "boolean" ? Number(v) : (v || "");
}

function render() {
  refreshChips();
  document.getElementById("gaps-only").setAttribute("aria-pressed", state.gapsOnly);
  document.getElementById("primary-only").setAttribute("aria-pressed", state.primaryOnly);
  document.getElementById("no-table-only").setAttribute("aria-pressed", state.noTableOnly);

  const rows = DATA.filter(matches).sort((a, b) => {
    const av = sortValue(a), bv = sortValue(b);
    return av < bv ? -state.sortDir : av > bv ? state.sortDir : 0;
  });

  document.querySelectorAll("thead th[data-key]").forEach(th => {
    th.removeAttribute("aria-sort");
    if (th.dataset.key === state.sortKey) th.setAttribute("aria-sort", state.sortDir === 1 ? "ascending" : "descending");
  });

  const tbody = document.getElementById("rows");
  tbody.innerHTML = rows.map(r => {
    const covered = Boolean(r.description);
    return "<tr>" +
      '<td class="col-cov"><span class="dot ' + (covered ? "good" : "gap") + '" title="' + (covered ? "Documented" : "Missing from metricDescriptions.yaml") + '"></span></td>' +
      '<td class="col-cov">' + (r.is_primary ? '<span class="star" title="Primary metric">★</span>' : "") + "</td>" +
      '<td class="col-cov">' + tableIcon(r.has_table) + "</td>" +
      '<td><div class="metric-name">' + escapeHtml(r.metric_name) + "</div></td>" +
      '<td><div class="task-name">' + escapeHtml(r.task_name) + '</div><div class="task-class">' + escapeHtml(r.task_class) +
        '</div><div class="task-atool">' + escapeHtml(r.atool) + "</div></td>" +
      "<td>" + bandDots(r.bands) + "</td>" +
      "<td>" + typeTags(r.metricTypes) + "</td>" +
      '<td class="units-sided">' + [escapeHtml(r.units), r.sided ? r.sided + "-sided" : ""].filter(Boolean).join(" · ") + "</td>" +
      '<td class="description' + (covered ? "" : " empty") + '">' + escapeHtml(r.description) + "</td>" +
      "</tr>";
  }).join("");

  document.getElementById("count-line").textContent = rows.length + " of " + DATA.length + " metrics shown";
  document.getElementById("empty-state").hidden = rows.length !== 0;
}

document.getElementById("search").addEventListener("input", e => { state.search = e.target.value.trim().toLowerCase(); render(); });
document.getElementById("gaps-only").addEventListener("click", () => { state.gapsOnly = !state.gapsOnly; render(); });
document.getElementById("primary-only").addEventListener("click", () => { state.primaryOnly = !state.primaryOnly; render(); });
document.getElementById("no-table-only").addEventListener("click", () => { state.noTableOnly = !state.noTableOnly; render(); });
document.querySelectorAll("thead th[data-key]").forEach(th => {
  th.addEventListener("click", () => {
    const key = th.dataset.key;
    if (state.sortKey === key) state.sortDir *= -1; else { state.sortKey = key; state.sortDir = 1; }
    render();
  });
});

render();
</script>
"""


def generate_html(df: pd.DataFrame, pipeline_path: str) -> str:
    """Render the metric table as a self-contained, filterable HTML page.

    Parameters
    ----------
    df
        Output of `build_metric_table`.
    pipeline_path
        Path to the drp_pipe pipeline YAML the table was built from, used to
        derive the page title and identity line.

    Returns
    -------
    str
        Full HTML document text.
    """
    instrument = Path(pipeline_path).parent.name
    pipeline_stem = Path(pipeline_path).stem
    pipeline_label = f"{instrument}/{Path(pipeline_path).name}"

    total = len(df)
    documented = int((df["description"] != "").sum())
    gaps = total - documented
    coverage_pct = round(100 * documented / total) if total else 0
    primary_count = int(df["is_primary"].sum())
    no_table_count = int((~df["has_table"]).sum())

    records = df.to_dict(orient="records")
    data_json = json.dumps(records, ensure_ascii=False).replace("</", "<\\/")
    band_colors_json = json.dumps(_BAND_COLORS)

    title = f"{instrument} {pipeline_stem} Metrics"
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    html = _HTML_TEMPLATE
    for placeholder, value in [
        ("__TITLE__", title),
        ("__H1__", f"{pipeline_stem} Metric Catalog"),
        ("__PIPELINE_LABEL__", pipeline_label),
        ("__TOTAL__", str(total)),
        ("__DOCUMENTED__", str(documented)),
        ("__GAPS__", str(gaps)),
        ("__COVERAGE_PCT__", str(coverage_pct)),
        ("__PRIMARY_COUNT__", str(primary_count)),
        ("__NO_TABLE__", str(no_table_count)),
        ("__GENERATED__", generated),
        ("__DATA_JSON__", data_json),
        ("__BAND_COLORS_JSON__", band_colors_json),
    ]:
        html = html.replace(placeholder, value)
    return html


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--pipeline",
        required=True,
        help="Path to the drp_pipe pipeline YAML, e.g. $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml.",
    )
    parser.add_argument(
        "--output",
        help="Write the metric table to this CSV path instead of printing it.",
    )
    parser.add_argument(
        "--html",
        help="Also write a self-contained, searchable/filterable HTML page to this path.",
    )
    parser.add_argument("--log-level", default="WARNING", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    log.info("Parsing pipeline: %s", args.pipeline)
    df = build_metric_table(args.pipeline)

    if df.empty:
        log.error("No metrics found in pipeline graph.")
        return 1

    if args.html:
        Path(args.html).write_text(generate_html(df, args.pipeline))
        log.info("Wrote HTML page to %s", args.html)

    if args.output:
        df.to_csv(args.output, index=False)
        log.info("Wrote %d metrics to %s", len(df), args.output)
    elif not args.html:
        with pd.option_context("display.max_rows", None, "display.max_colwidth", 60):
            print(df.to_string(index=False))
        print(f"\n{len(df)} metrics.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
