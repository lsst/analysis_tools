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
from __future__ import annotations

from unittest import TestCase, main

import numpy as np

import lsst.utils.tests
from lsst.analysis.tools.actions.scalar import CountAction, MaxAction, MedianAction, MinAction
from lsst.analysis.tools.atools import AggregatedTaskMetadataMetricTool
from lsst.analysis.tools.tasks import AggregatedTaskMetadataAnalysisTask
from lsst.analysis.tools.tasks.metadataAnalysis import AggregatedTaskMetadataAnalysisConfig
from lsst.pipe.base import UpstreamFailureNoWorkFound

_PLOT_INFO = {"tableName": "", "run": "", "plotName": ""}


def _make_tool(metrics, aggregation_units=None, subtask_names=None):
    """Create and finalize a minimal AggregatedTaskMetadataMetricTool.

    Parameters
    ----------
    metrics : `dict` [`str`, `str`]
        Metric names and their units.
    aggregation_units : `dict` [`str`, `str`], optional
        Fixed units keyed by aggregation suffix.
    subtask_names : `dict` [`str`, `str`], optional
        Subtask names keyed by metric name.

    Returns
    -------
    tool : `AggregatedTaskMetadataMetricTool`
        A finalized tool.
    """
    tool = AggregatedTaskMetadataMetricTool()
    tool.taskName = "myTask"
    tool.metrics = metrics
    if subtask_names is not None:
        tool.subTaskNames = subtask_names
    if aggregation_units is not None:
        tool.aggregationUnits = aggregation_units

    aggregations = {"min": MinAction, "max": MaxAction, "median": MedianAction, "ct": CountAction}
    for metric in tool.metrics:
        sanitized = metric.replace(" ", "_")
        for agg_name, agg_cls in aggregations.items():
            action = agg_cls()
            action.vectorKey = metric
            setattr(tool.process.calculateActions, f"{sanitized}_{agg_name}", action)

    tool.finalize()
    return tool


def _make_task(metrics, subtask_names=None):
    """Create a minimal AggregatedTaskMetadataAnalysisTask.

    Parameters
    ----------
    metrics : `dict` [`str`, `str`]
        Metric names and their units.
    subtask_names : `dict` [`str`, `str`], optional
        Subtask names keyed by metric name.

    Returns
    -------
    task : `AggregatedTaskMetadataAnalysisTask`
        A configured task.
    """
    config = AggregatedTaskMetadataAnalysisConfig()
    config.connections.inputName = "myTask_metadata"
    config.connections.outputName = "myTask_metadata_agg"
    config.atools.myTool = AggregatedTaskMetadataMetricTool
    config.atools.myTool.taskName = "myTask"
    config.atools.myTool.metrics = metrics
    if subtask_names is not None:
        config.atools.myTool.subTaskNames = subtask_names

    action = MedianAction()
    action.vectorKey = next(iter(metrics))
    config.atools.myTool.process.calculateActions.metric_median = action

    return AggregatedTaskMetadataAnalysisTask(config=config)


class _MockHandle:
    """Minimal stand-in for a deferred dataset handle."""

    def __init__(self, metadata_dict):
        self._metadata = metadata_dict

    def get(self):
        metadata = self._metadata

        class _MockMetadata:
            def to_dict(self):
                return metadata

        return _MockMetadata()


def _make_handles(metadata_dicts):
    """Create a list of mock handles from a list of metadata dicts.

    Parameters
    ----------
    metadata_dicts : `list` [`dict`]
        One dict per handle, as returned by ``TaskMetadata.to_dict()``.

    Returns
    -------
    handles : `list` [`_MockHandle`]
    """
    return [_MockHandle(md) for md in metadata_dicts]


class TestAggregatedTaskMetadataMetricToolFinalize(TestCase):
    """Tests for AggregatedTaskMetadataMetricTool.finalize."""

    def testUnitsInheritedFromSourceMetric(self):
        """min/max/median aggregations should inherit the source metric's
        unit."""
        tool = _make_tool({"nStars": "ct", "psfSigma": "pixel"})
        units = tool.produce.metric.units
        self.assertEqual(units["nStars_min"], "ct")
        self.assertEqual(units["nStars_max"], "ct")
        self.assertEqual(units["nStars_median"], "ct")
        self.assertEqual(units["psfSigma_min"], "pixel")
        self.assertEqual(units["psfSigma_max"], "pixel")
        self.assertEqual(units["psfSigma_median"], "pixel")

    def testAggregationUnitsOverride(self):
        """aggregationUnits should override the inherited unit for that
        suffix."""
        tool = _make_tool(
            {"nStars": "ct", "psfSigma": "pixel"},
            aggregation_units={"ct": "ct"},
        )
        units = tool.produce.metric.units
        # Count aggregations should use the override unit regardless of source.
        self.assertEqual(units["nStars_ct"], "ct")
        self.assertEqual(units["psfSigma_ct"], "ct")
        # Other aggregations should still inherit.
        self.assertEqual(units["psfSigma_min"], "pixel")

    def testSpacesInMetricNames(self):
        """Metric names containing spaces must be handled gracefully."""
        tool = _make_tool({"n stars": "ct"})
        units = tool.produce.metric.units
        # The action name uses underscores; units should still be populated.
        self.assertIn("n_stars_min", units)
        self.assertIn("n_stars_max", units)
        self.assertIn("n_stars_median", units)
        self.assertEqual(units["n_stars_min"], "ct")

    def testActionWithUnknownVectorKeyRaisesError(self):
        """An action whose vectorKey is not in metrics should raise ValueError.

        This guards against misconfiguration where someone assigns a vectorKey
        that doesn't correspond to any configured metric.
        """
        tool = AggregatedTaskMetadataMetricTool()
        tool.taskName = "myTask"
        tool.metrics = {"nStars": "ct"}

        known_action = MedianAction()
        known_action.vectorKey = "nStars"
        tool.process.calculateActions.nStars_median = known_action

        unknown_action = MedianAction()
        unknown_action.vectorKey = "notAMetric"
        tool.process.calculateActions.notAMetric_median = unknown_action

        with self.assertRaises(ValueError):
            tool.finalize()


class TestAggregatedTaskMetadataMetricToolNumerics(TestCase):
    """Numerical correctness of aggregations."""

    def setUp(self):
        self.values = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        self.tool = _make_tool({"nStars": "ct"}, aggregation_units={"ct": "ct"})

    def _call_tool(self, data):
        return self.tool(data, plotInfo=dict(_PLOT_INFO))

    def testMin(self):
        result = self._call_tool({"nStars": self.values})
        self.assertAlmostEqual(result["nStars_min"].quantity.value, 10.0)

    def testMax(self):
        result = self._call_tool({"nStars": self.values})
        self.assertAlmostEqual(result["nStars_max"].quantity.value, 50.0)

    def testMedian(self):
        result = self._call_tool({"nStars": self.values})
        self.assertAlmostEqual(result["nStars_median"].quantity.value, 30.0)

    def testCount(self):
        result = self._call_tool({"nStars": self.values})
        self.assertAlmostEqual(result["nStars_ct"].quantity.value, 5.0)

    def testCountExcludesNaN(self):
        """Count should exclude NaN values."""
        values_with_nan = np.array([10.0, np.nan, 30.0, np.nan, 50.0])
        result = self._call_tool({"nStars": values_with_nan})
        self.assertAlmostEqual(result["nStars_ct"].quantity.value, 3.0)


class TestAggregatedTaskMetadataAnalysisTask(TestCase):
    """Tests for AggregatedTaskMetadataAnalysisTask._collectData."""

    def testEmptyMetadataRaisesNoWorkFound(self):
        """All inputs returning empty metadata should raise NoWorkFound."""
        task = _make_task({"nStars": "ct"})
        handles = _make_handles([{}, {}, {}])
        with self.assertRaises(UpstreamFailureNoWorkFound):
            task._collectData(handles, "myTask")

    def testMissingMetricRaisesNoWorkFound(self):
        """A configured metric absent from all inputs should raise
        NoWorkFound."""
        task = _make_task({"nStars": "ct", "psfSigma": "pixel"})
        # Inputs only contain nStars, not psfSigma.
        handles = _make_handles([
            {"myTask": {"nStars": 10.0}},
            {"myTask": {"nStars": 20.0}},
        ])
        with self.assertRaises(UpstreamFailureNoWorkFound):
            task._collectData(handles, "myTask")

    def testPartialInputsAccepted(self):
        """Inputs missing a metric in some (but not all) detectors are fine."""
        task = _make_task({"nStars": "ct"})
        # One input lacks nStars — it should be silently skipped.
        handles = _make_handles([
            {"myTask": {"nStars": 10.0}},
            {"myTask": {}},
            {"myTask": {"nStars": 30.0}},
        ])
        data = task._collectData(handles, "myTask")
        self.assertEqual(len(data["nStars"]), 2)
        self.assertIn(10.0, data["nStars"])
        self.assertIn(30.0, data["nStars"])
        # Verify values are collected as-is, not transformed.
        self.assertAlmostEqual(min(data["nStars"]), 10.0)
        self.assertAlmostEqual(max(data["nStars"]), 30.0)

    def testSubTaskNamesRespected(self):
        """Metrics in subtasks should be found under the correct key."""
        task = _make_task(
            {"cosmicRayCount": "ct"},
            subtask_names={"cosmicRayCount": "repair"},
        )
        handles = _make_handles([
            {"myTask:repair": {"cosmicRayCount": 5.0}},
            {"myTask:repair": {"cosmicRayCount": 8.0}},
        ])
        data = task._collectData(handles, "myTask")
        self.assertIn("cosmicRayCount", data)
        self.assertEqual(len(data["cosmicRayCount"]), 2)

    def testTypoInMetricNameRaisesNoWorkFound(self):
        """A metric name that doesn't match any key in the metadata raises
        NoWorkFound, covering the case of a configuration typo.
        """
        task = _make_task({"nSatrs": "ct"})  # typo: nSatrs instead of nStars
        handles = _make_handles([
            {"myTask": {"nStars": 10.0}},
            {"myTask": {"nStars": 20.0}},
        ])
        with self.assertRaises(UpstreamFailureNoWorkFound):
            task._collectData(handles, "myTask")


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
