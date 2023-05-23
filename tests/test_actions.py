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

import unittest

import astropy.units as u
import numpy as np
import pandas as pd
from lsst.analysis.tools.actions.scalar.scalarActions import (
    ApproxFloor,
    CountAction,
    MeanAction,
    MedianAction,
    SigmaMadAction,
    StdevAction,
)
from lsst.analysis.tools.actions.vector.calcBinnedStats import CalcBinnedStatsAction
from lsst.analysis.tools.actions.vector.calcShapeSize import CalcShapeSize
from lsst.analysis.tools.actions.vector.mathActions import (
    AddVector,
    ConstantValue,
    DivideVector,
    FractionalDifference,
    MultiplyVector,
    SubtractVector,
)
from lsst.analysis.tools.actions.vector.selectors import (
    CoaddPlotFlagSelector,
    FlagSelector,
    GalaxySelector,
    RangeSelector,
    SkyObjectSelector,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from lsst.analysis.tools.actions.vector.vectorActions import (
    ConvertFluxToMag,
    DownselectVector,
    ExtinctionCorrectedMagDiff,
    LoadVector,
    MagDiff,
)


class TestScalarActions(unittest.TestCase):
    """ "Test ScalarActions"""

    def setUp(self):
        x = np.arange(100, dtype=float)
        x[31] = np.nan
        x[41] = np.nan
        x[59] = np.nan
        y = x**2
        self.data = {"r_y": x, "i_y": y}
        self.mask = np.zeros(100, dtype=bool)
        self.mask[1:50] = 1

    def _testScalarActionAlmostEqual(self, cls, truth, maskedTruth):
        action = cls(vectorKey="{band}_y")
        schema = [col for col, colType in action.getInputSchema()]
        self.assertEqual(schema, ["{band}_y"])
        result = action(self.data, band="i")
        self.assertAlmostEqual(result, truth)
        result = action(self.data, mask=self.mask, band="i")
        self.assertAlmostEqual(result, maskedTruth)

    def testMedianAction(self):
        self._testScalarActionAlmostEqual(MedianAction, 2500, 576)

    def testMeanAction(self):
        self._testScalarActionAlmostEqual(MeanAction, 3321.9278350515465, 803.8936170212766)

    def testStdevAction(self):
        self._testScalarActionAlmostEqual(StdevAction, 2984.5855649976297, 733.5989754407842)

    def testSigmaMadAction(self):
        self._testScalarActionAlmostEqual(SigmaMadAction, 3278.033505115886, 759.0923358748682)

    def testCountAction(self):
        self._testScalarActionAlmostEqual(CountAction, 97, 47)

    def testApproxFloorAction(self):
        self._testScalarActionAlmostEqual(ApproxFloor, 9216.0, 2352.5)


class TestVectorActions(unittest.TestCase):
    """Test VectorActions"""

    def setUp(self):
        self.size = 5
        r = np.arange(self.size) + 1
        i = r**2
        rFlag = np.ones(self.size)
        rFlag[1] = 0
        iFlag = np.ones(self.size)
        iFlag[3] = 0
        data = {
            "r_vector": r,
            "i_vector": i,
            "r_flag": rFlag,
            "i_flag": iFlag,
            "g_vector": 3 * r,
            "E(B-V)": r[::-1],
            "r_ixx": r**2,
            "r_iyy": r[::-1] ** 2,
            "r_ixy": r * r[::-1],
            "g_iixx": r**2,
            "g_iiyy": r[::-1] ** 2,
            "g_iixy": r * r[::-1],
        }

        self.data = pd.DataFrame.from_dict(data)

    def _checkSchema(self, action, truth):
        schema = sorted([col for col, colType in action.getInputSchema()])
        self.assertEqual(schema, truth)

    # VectorActions with their own files

    def testCalcBinnedStats(self):
        selector = RangeSelector(key="r_vector", minimum=0, maximum=self.size + 1)
        prefix = "a_"
        stats = CalcBinnedStatsAction(name_prefix=prefix, selector_range=selector, key_vector="r_vector")
        result = stats(self.data)
        median = (1 + self.size) / 2.0
        truth = {
            stats.name_mask: np.ones(self.size),
            stats.name_median: median,
            stats.name_sigmaMad: 1.482602218505602 * np.median(np.abs(self.data["r_vector"] - median)),
            stats.name_count: self.size,
            stats.name_select_maximum: self.size,
            stats.name_select_median: median,
            stats.name_select_minimum: 1,
            "range_maximum": self.size + 1,
            "range_minimum": 0,
        }
        self.assertEqual(list(result.keys()), list(truth.keys()))

        self.assertAlmostEqual(result[stats.name_sigmaMad], truth[stats.name_sigmaMad])
        del truth[stats.name_sigmaMad]

        np.testing.assert_array_equal(result[stats.name_mask], truth[stats.name_mask])
        del truth[stats.name_mask]

        for key, value in truth.items():
            self.assertEqual(result[key], value, key)

    # def testCalcRhoStatistics(self): TODO: implement

    def testCalcShapeSize(self):
        xx = self.data["r_ixx"]
        yy = self.data["r_iyy"]
        xy = self.data["r_ixy"]

        # Test determinant with defaults
        action = CalcShapeSize()
        result = action(self.data, band="r")
        schema = [col for col, colType in action.getInputSchema()]
        self.assertEqual(sorted(schema), ["{band}_ixx", "{band}_ixy", "{band}_iyy"])
        truth = 0.25 * (xx * yy - xy**2)
        np.testing.assert_array_almost_equal(result, truth)

        # Test trace with columns specified
        action = CalcShapeSize(
            colXx="{band}_iixx",
            colYy="{band}_iiyy",
            colXy="{band}_iixy",
            sizeType="trace",
        )
        result = action(self.data, band="g")
        schema = [col for col, colType in action.getInputSchema()]
        self.assertEqual(sorted(schema), ["{band}_iixx", "{band}_iiyy"])
        truth = np.sqrt(0.5 * (xx + yy))
        np.testing.assert_array_almost_equal(result, truth)

    # def testCalcE(self): TODO: implement

    # def testCalcEDiff(self): TODO: implement

    # def testCalcE1(self): TODO: implement

    # def testCalcE2(self): TODO: implement

    # MathActions

    def _testMath(self, ActionType, truth, compare_exact: bool = False):
        actionA = LoadVector(vectorKey="{band1}_vector")
        actionB = LoadVector(vectorKey="{band2}_vector")
        action = ActionType(actionA=actionA, actionB=actionB)
        result = action(self.data, band1="r", band2="i")
        self._checkSchema(action, [actionA.vectorKey, actionB.vectorKey])
        if compare_exact:
            np.testing.assert_array_equal(result, truth)
        else:
            np.testing.assert_array_almost_equal(result, truth)

    def testConstant(self):
        truth = [42.0]
        action = ConstantValue(value=truth[0])
        self._checkSchema(action, [])
        result = action({})
        np.testing.assert_array_equal(result, truth)

    def testAdd(self):
        truth = [2.0, 6.0, 12.0, 20.0, 30.0]
        self._testMath(AddVector, truth, True)

    def testSubtract(self):
        truth = [0.0, -2.0, -6.0, -12.0, -20.0]
        self._testMath(SubtractVector, truth, True)

    def testMultiply(self):
        truth = [1.0, 8.0, 27.0, 64.0, 125.0]
        self._testMath(MultiplyVector, truth, False)

    def testDivide(self):
        truth = 1 / np.arange(1, 6)
        self._testMath(DivideVector, truth, False)

    def testFractionalDifference(self):
        actionA = LoadVector(vectorKey="{band1}_vector")
        actionB = LoadVector(vectorKey="{band2}_vector")
        truth = [0.0, -0.5, -0.6666666666666666, -0.75, -0.8]
        diff = FractionalDifference(actionA=actionA, actionB=actionB)
        result = diff(self.data, band1="r", band2="i")
        self._checkSchema(diff, ["{band1}_vector", "{band2}_vector"])
        np.testing.assert_array_almost_equal(result, truth)

    # Basic vectorActions

    def testLoadVector(self):
        action = LoadVector(vectorKey="vector")
        inputData = {"vector": np.arange(5)}
        result = action(inputData)
        np.testing.assert_array_equal(result, inputData["vector"])

    def testDownselectVector(self):
        selector = FlagSelector(selectWhenTrue=["{band}_flag"])
        action = DownselectVector(vectorKey="{band}_vector", selector=selector)
        result = action(self.data, band="r")
        self._checkSchema(action, ["{band}_flag", "{band}_vector"])
        np.testing.assert_array_equal(result, np.array([1, 3, 4, 5]))

    # def testMultiCriteriaDownselectVector(self): TODO: implement

    # Astronomical vectorActions

    # def testCalcSn(self): TODO: implement

    def testConvertFluxToMag(self):
        truth = [
            31.4,
            29.89485002168,
            29.0143937264,
            28.38970004336,
            27.90514997832,
        ]
        action = ConvertFluxToMag(vectorKey="{band}_vector")
        result = action(self.data, band="i")
        self._checkSchema(action, ["{band}_vector"])
        np.testing.assert_array_almost_equal(result, truth)

    # def testConvertUnits(self): TODO: implement

    def testMagDiff(self):
        # Use the same units as the data so that we know the difference exactly
        # without conversions.
        # testExtinctionCorrectedMagDiff will test the conversions
        truth = np.array(2 * self.data["r_vector"]) * u.ABmag
        action = MagDiff(
            col1="{band1}_vector",
            col2="{band2}_vector",
            fluxUnits1="mag(AB)",
            fluxUnits2="mag(AB)",
            returnMillimags=False,
        )
        result = action(self.data, band1="g", band2="r")
        self._checkSchema(action, ["{band1}_vector", "{band2}_vector"])
        np.testing.assert_array_almost_equal(result, truth.value)

    def testExtinctionCorrectedMagDiff(self):
        for returnMillimags in (True, False):
            # Check that conversions are working properly
            magDiff = MagDiff(
                col1="{band1}_vector",
                col2="{band2}_vector",
                fluxUnits1="jansky",
                fluxUnits2="jansky",
                returnMillimags=returnMillimags,
            )
            action = ExtinctionCorrectedMagDiff(
                magDiff=magDiff,
                band1="g",
                band2="r",
                ebvCol="E(B-V)",
                extinctionCoeffs={"g": 0.2, "r": 1.5},
            )

            result = action(self.data, band1="g", band2="r")
            lhs = (np.array(self.data["g_vector"]) * u.jansky).to(u.ABmag)
            rhs = (np.array(self.data["r_vector"]) * u.jansky).to(u.ABmag)
            diff = lhs - rhs
            correction = np.array((0.2 - 1.5) * self.data["E(B-V)"]) * u.mag
            if returnMillimags:
                diff = diff.to(u.mmag)
                correction = correction.to(u.mmag)
            truth = diff - correction
            self._checkSchema(action, ["E(B-V)", "{band1}_vector", "{band2}_vector"])
            np.testing.assert_array_almost_equal(result, truth.value)

        # Test with hard coded bands
        magDiff = MagDiff(
            col1="g_vector",
            col2="r_vector",
            fluxUnits1="jansky",
            fluxUnits2="jansky",
            returnMillimags=False,
        )
        action = ExtinctionCorrectedMagDiff(
            magDiff=magDiff,
            ebvCol="E(B-V)",
            extinctionCoeffs={"g": 0.2, "r": 1.5},
        )
        result = action(self.data)
        lhs = (np.array(self.data["g_vector"]) * u.jansky).to(u.ABmag)
        rhs = (np.array(self.data["r_vector"]) * u.jansky).to(u.ABmag)
        diff = lhs - rhs
        correction = np.array((0.2 - 1.5) * self.data["E(B-V)"]) * u.mag
        truth = diff - correction
        self._checkSchema(action, ["E(B-V)", "g_vector", "r_vector"])
        np.testing.assert_array_almost_equal(result, truth.value)

    # def testRAcosDec(self): TODO: implement

    # Statistical vectorActions

    # def testPerGroupStatistic(self): TODO: implement

    # def testResidualWithPerGroupStatistic(self): TODO: implement


class TestVectorSelectors(unittest.TestCase):
    def setUp(self):
        self.size = 20
        falseFlags = {
            "{band}_psfFlux_flag": [1],
            "{band}_pixelFlags_saturatedCenter": [3],
            "{band}_extendedness_flag": [5],
            "xy_flag": [7],
            "i_pixelFlags_edge": [13],
            "r_pixelFlags_edge": [15],
            "sky_object": [13, 15, 17],
        }

        trueFlags = {
            "detect_isPatchInner": [9],
            "detect_isDeblendedSource": [11],
        }

        flux = np.arange(self.size) * 10
        fluxErr = np.ones(self.size) * 0.1
        extendedness = np.arange(20) / 20 - 0.1

        self.data = {
            "r_psfFlux": flux,
            "r_psfFluxErr": fluxErr,
            "i_cmodelFlux": flux[::-1],
            "i_cmodelFluxError": fluxErr[::-1],
            "r_cmodelFlux": flux,
            "r_cmodelFluxError": fluxErr,
            "i_extendedness": extendedness,
            "i_extended": extendedness,
        }
        bands = ("r", "i")
        for band in bands:
            for flag, bits in falseFlags.items():
                vector = np.zeros(self.size, dtype=bool)
                for bit in bits:
                    vector[bit] = 1
                self.data[flag.format(band=band)] = vector

            for flag, bits in trueFlags.items():
                vector = np.ones(self.size, dtype=bool)
                for bit in bits:
                    vector[bit] = 0
                self.data[flag.format(band=band)] = vector

    def _checkSchema(self, action, truth):
        schema = [col for col, colType in action.getInputSchema()]
        self.assertEqual(sorted(schema), sorted(truth))

    def testFlagSelector(self):
        selector = FlagSelector(
            selectWhenFalse=["{band}_psfFlux_flag"], selectWhenTrue=["detect_isPatchInner"]
        )
        self._checkSchema(selector, ["detect_isPatchInner", "{band}_psfFlux_flag"])
        result = selector(self.data, band="r")
        truth = np.ones(self.size, dtype=bool)
        truth[1] = False
        truth[9] = False
        np.testing.assert_array_equal(result, truth)

    def testCoaddPlotFlagSelector(self):
        # Test defaults
        # Bands needs to be set to something otherwise it
        # will crash as the default value looks it up
        # elsewhere.
        selector = CoaddPlotFlagSelector(bands=["i"])
        keys = [
            "{band}_psfFlux_flag",
            "{band}_pixelFlags_saturatedCenter",
            "{band}_extendedness_flag",
            "xy_flag",
            "detect_isPatchInner",
            "detect_isDeblendedSource",
        ]
        self._checkSchema(selector, keys)

        result = selector(self.data)
        truth = np.ones(self.size, dtype=bool)
        for bit in (1, 3, 5, 7, 9, 11):
            truth[bit] = 0
        np.testing.assert_array_equal(result, truth)

        # Test bands override
        selector = CoaddPlotFlagSelector(
            bands=["i", "r"],
            selectWhenFalse=["{band}_psfFlux_flag"],
            selectWhenTrue=["detect_isDeblendedSource"],
        )
        self._checkSchema(selector, ["{band}_psfFlux_flag", "detect_isDeblendedSource"])
        result = selector(self.data)
        truth = np.ones(self.size, dtype=bool)
        for bit in (1, 11):
            truth[bit] = 0
        np.testing.assert_array_equal(result, truth)

    def testRangeSelector(self):
        selector = RangeSelector(key="r_psfFlux", minimum=np.nextafter(20, 30), maximum=50)
        self._checkSchema(selector, ["r_psfFlux"])
        result = self.data["r_psfFlux"][selector(self.data)]
        truth = [30, 40]
        np.testing.assert_array_equal(result, truth)

    def testSnSelector(self):
        # test defaults
        selector = SnSelector()
        keys = [
            "{band}_psfFlux",
            "{band}_psfFluxErr",
        ]
        self._checkSchema(selector, keys)
        result = selector(self.data, bands=["r"])
        truth = np.ones(self.size, dtype=bool)
        truth[:6] = 0
        np.testing.assert_array_equal(result, truth)

        # test overrides
        selector = SnSelector(
            fluxType="{band}_cmodelFlux",
            threshold=200.0,
            uncertaintySuffix="Error",
            bands=["r", "i"],
        )
        keys = [
            "{band}_cmodelFlux",
            "{band}_cmodelFluxError",
        ]
        self._checkSchema(selector, keys)
        result = selector(self.data)
        truth = np.ones(self.size, dtype=bool)
        truth[:3] = 0
        truth[-3:] = 0
        np.testing.assert_array_equal(result, truth)

    def testSkyObjectSelector(self):
        # Test default
        selector = SkyObjectSelector()
        keys = ["{band}_pixelFlags_edge", "sky_object"]
        self._checkSchema(selector, keys)
        result = selector(self.data)
        truth = np.zeros(self.size, dtype=bool)
        truth[15] = 1
        truth[17] = 1
        np.testing.assert_array_equal(result, truth)

        # Test overrides
        selector = SkyObjectSelector(bands=["i", "r"])
        self._checkSchema(selector, keys)
        result = selector(self.data)
        truth = np.zeros(self.size, dtype=bool)
        truth[17] = 1
        np.testing.assert_array_equal(result, truth)

    def testStarSelector(self):
        # test default
        selector = StarSelector()
        self._checkSchema(selector, ["{band}_extendedness"])
        result = selector(self.data, band="i")
        truth = (self.data["i_extendedness"] >= 0) & (self.data["i_extendedness"] < 0.5)
        np.testing.assert_array_almost_equal(result, truth)

        # Test overrides
        selector = StarSelector(vectorKey="i_extended", extendedness_maximum=0.3)
        result = selector(self.data, band="i")
        truth = (self.data["i_extendedness"] >= 0) & (self.data["i_extendedness"] < 0.3)
        np.testing.assert_array_almost_equal(result, truth)

    def testGalaxySelector(self):
        # test default
        selector = GalaxySelector()
        self._checkSchema(selector, ["{band}_extendedness"])
        result = selector(self.data, band="i")
        truth = self.data["i_extendedness"] > 0.5
        np.testing.assert_array_almost_equal(result, truth)

        # Test overrides
        selector = GalaxySelector(vectorKey="i_extended", extendedness_minimum=0.3)
        result = selector(self.data, band="i")
        truth = self.data["i_extendedness"] > 0.3
        np.testing.assert_array_almost_equal(result, truth)

    def testVectorSelector(self):
        selector = VectorSelector(vectorKey="{band}_psfFlux_flag")
        self._checkSchema(selector, ["{band}_psfFlux_flag"])
        result = selector(self.data, band="i")
        truth = np.zeros(self.size, dtype=bool)
        truth[1] = True
        np.testing.assert_array_equal(result, truth)
