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
from lsst.analysis.tools.actions.scalar import (
    CountAction,
    MaxAction,
    MeanAction,
    MinAction,
    WeightedMeanAction,
)
from lsst.analysis.tools.atools import WholeTractMaskFractionTool
from lsst.analysis.tools.tasks import WholeTractMaskFractionAnalysisTask
from lsst.pex.exceptions import InvalidParameterError
from lsst.pipe.base import InMemoryDatasetHandle

_PLANE_BITS = {"NO_DATA": 1, "SAT": 2, "CR": 4}


class _MockMask:
    def __init__(self, array, plane_bits=_PLANE_BITS):
        self._array = array
        self._plane_bits = plane_bits

    @property
    def array(self):
        return self._array

    def getPlaneBitMask(self, plane):
        if plane not in self._plane_bits:
            raise InvalidParameterError(f"Unknown mask plane: {plane!r}")
        return self._plane_bits[plane]


class _MockExposure:
    def __init__(self, array, plane_bits=_PLANE_BITS):
        self._mask = _MockMask(array, plane_bits)

    def getMask(self):
        return self._mask


def _make_tool(mask_planes):
    """Create and finalize a WholeTractMaskFractionTool with standard
    aggregations."""
    tool = WholeTractMaskFractionTool()

    all_planes = set(mask_planes) | {"NO_DATA"}
    aggregations = {"min": MinAction, "max": MaxAction, "mean": MeanAction, "ct": CountAction}

    for plane in all_planes:
        for agg_name, agg_cls in aggregations.items():
            action = agg_cls()
            action.vectorKey = f"{plane}_fraction"
            setattr(tool.process.calculateActions, f"{agg_name}_{plane}_fraction", action)

    for plane in mask_planes:
        if plane != "NO_DATA":
            action = WeightedMeanAction()
            action.vectorKey = f"{plane}_valid_data_fraction"
            action.weightsKey = "valid_data_pixel_count"
            setattr(tool.process.calculateActions, f"mean_{plane}_valid_data_fraction", action)

    tool.finalize()
    return tool


def _make_task(mask_planes):
    """Create a WholeTractMaskFractionAnalysisTask with one configured tool."""
    config = WholeTractMaskFractionAnalysisTask.ConfigClass()
    config.connections.outputName = "wholeTractMaskFraction"
    config.maskPlanes = mask_planes
    config.atools.maskFractionMetrics = WholeTractMaskFractionTool
    return WholeTractMaskFractionAnalysisTask(config=config)


class TestWholeTractMaskFractionToolFinalize(TestCase):

    def testUnitsAreEmptyString(self):
        tool = _make_tool(["SAT", "CR"])
        for unit in tool.produce.metric.units.values():
            self.assertEqual(unit, "")


class TestWholeTractMaskFractionTaskCollectAllPlanes(TestCase):

    def testConfiguredPlanesIncluded(self):
        task = _make_task(["SAT", "CR"])
        self.assertIn("SAT", task.config.maskPlanes)
        self.assertIn("CR", task.config.maskPlanes)


class TestWholeTractMaskFractionTaskComputeKeyedData(TestCase):

    def setUp(self):
        self.task = _make_task(["SAT", "CR"])

    def _compute(self, handles, planes=None):
        if planes is None:
            planes = {"SAT", "CR", "NO_DATA"}
        return self.task._computeMaskFractions(handles, planes)

    def testFractionAllFlagged(self):
        """All pixels flagged → fraction == 1.0."""
        arr = np.full((4, 4), _PLANE_BITS["SAT"], dtype=np.int32)
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertAlmostEqual(result["SAT_fraction"][0], 1.0)

    def testFractionNoneFlagged(self):
        """No pixels flagged → fraction == 0.0."""
        arr = np.zeros((4, 4), dtype=np.int32)
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertAlmostEqual(result["SAT_fraction"][0], 0.0)

    def testFractionHalfFlagged(self):
        arr = np.zeros((4, 4), dtype=np.int32)
        arr[:2, :] = _PLANE_BITS["SAT"]  # half the pixels
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertAlmostEqual(result["SAT_fraction"][0], 0.5)

    def testNoDataAlwaysIncluded(self):
        arr = np.zeros((4, 4), dtype=np.int32)
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertIn("NO_DATA_fraction", result)

    def testNoDataHasNoFractionValid(self):
        arr = np.zeros((4, 4), dtype=np.int32)
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertNotIn("NO_DATA_valid_data_fraction", result)

    def testFractionValidExcludesNoDataPixels(self):
        # 8 pixels total: 4 NO_DATA, 4 valid; of the 4 valid, 2 are SAT
        arr = np.zeros((2, 4), dtype=np.int32)
        arr[0, :] = _PLANE_BITS["NO_DATA"]  # row 0: NO_DATA
        arr[1, :2] = _PLANE_BITS["SAT"]  # row 1: 2 SAT (valid)
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertAlmostEqual(result["SAT_valid_data_fraction"][0], 0.5)

    def testValidPixelCountCorrect(self):
        arr = np.zeros((2, 4), dtype=np.int32)
        arr[0, :] = _PLANE_BITS["NO_DATA"]  # 4 NO_DATA pixels
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertEqual(result["valid_data_pixel_count"][0], 4)

    def testUnknownPlaneSkipped(self):
        arr = np.zeros((4, 4), dtype=np.int32)
        planes = {"SAT", "CR", "NO_DATA", "UNKNOWN_PLANE"}
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))], planes=planes)
        self.assertNotIn("UNKNOWN_PLANE_fraction", result)

    def testMultiplePatchesAccumulate(self):
        arr_a = np.zeros((4, 4), dtype=np.int32)
        arr_b = np.full((4, 4), _PLANE_BITS["SAT"], dtype=np.int32)
        result = self._compute(
            [InMemoryDatasetHandle(_MockExposure(arr_a)), InMemoryDatasetHandle(_MockExposure(arr_b))]
        )
        self.assertEqual(len(result["SAT_fraction"]), 2)
        self.assertAlmostEqual(result["SAT_fraction"][0], 0.0)
        self.assertAlmostEqual(result["SAT_fraction"][1], 1.0)

    def testValidPixelCountAcrossMultiplePatches(self):
        arr_a = np.zeros((2, 4), dtype=np.int32)
        arr_a[0, :] = _PLANE_BITS["NO_DATA"]  # 4 valid pixels
        arr_b = np.zeros((2, 4), dtype=np.int32)  # 8 valid pixels
        result = self._compute(
            [InMemoryDatasetHandle(_MockExposure(arr_a)), InMemoryDatasetHandle(_MockExposure(arr_b))]
        )
        np.testing.assert_array_equal(result["valid_data_pixel_count"], [4, 8])

    def testAllNoDataPatchYieldsNaNFractionValid(self):
        """A patch where all pixels are NO_DATA should yield NaN for
        _fraction_valid, not be omitted, so the vector length matches
        _fraction."""
        arr_all_no_data = np.full((4, 4), _PLANE_BITS["NO_DATA"], dtype=np.int32)
        arr_normal = np.zeros((4, 4), dtype=np.int32)
        arr_normal[:2, :] = _PLANE_BITS["SAT"]
        result = self._compute(
            [
                InMemoryDatasetHandle(_MockExposure(arr_all_no_data)),
                InMemoryDatasetHandle(_MockExposure(arr_normal)),
            ]
        )
        # Both patches contribute to both vectors
        self.assertEqual(len(result["SAT_fraction"]), 2)
        self.assertEqual(len(result["SAT_valid_data_fraction"]), 2)
        # All-NO_DATA patch → NaN; normal patch → 0.5
        self.assertTrue(np.isnan(result["SAT_valid_data_fraction"][0]))
        self.assertAlmostEqual(result["SAT_valid_data_fraction"][1], 0.5)

    def testTwoPlanesComputedIndependently(self):
        """SAT and CR fractions should be computed independently."""
        arr = np.zeros((4, 4), dtype=np.int32)
        arr[:1, :] = _PLANE_BITS["SAT"]  # 4/16 = 0.25 SAT
        arr[1:2, :] = _PLANE_BITS["CR"]  # 4/16 = 0.25 CR
        result = self._compute([InMemoryDatasetHandle(_MockExposure(arr))])
        self.assertAlmostEqual(result["SAT_fraction"][0], 0.25)
        self.assertAlmostEqual(result["CR_fraction"][0], 0.25)


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
