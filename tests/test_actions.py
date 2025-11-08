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
import lsst.utils.tests
import numpy as np
import pandas as pd
from lsst.analysis.tools.actions.keyedData.calcDistances import CalcRelativeDistances
from lsst.analysis.tools.actions.scalar import (
    ApproxFloor,
    CountAction,
    MeanAction,
    MedianAction,
    SigmaMadAction,
    StdevAction,
)
from lsst.analysis.tools.actions.vector import (
    CalcBinnedStatsAction,
    CalcMomentSize,
    CalcRhoStatistics,
    ConvertFluxToMag,
    DownselectVector,
    ExtinctionCorrectedMagDiff,
    LoadVector,
    MagDiff,
)
from lsst.analysis.tools.actions.vector.mathActions import (
    AddVector,
    ConstantValue,
    DivideVector,
    FractionalDifference,
    Log10Vector,
    MultiplyVector,
    RaiseFromBaseVector,
    RaiseToPowerVector,
    SqrtVector,
    SquareVector,
    SubtractVector,
)
from lsst.analysis.tools.actions.vector.selectors import (
    CoaddPlotFlagSelector,
    FlagSelector,
    GalaxySelector,
    RangeSelector,
    SetSelector,
    SkyObjectSelector,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from lsst.analysis.tools.math import sqrt
from lsst.pex.config import FieldValidationError


class TestScalarActions(unittest.TestCase):
    """ "Test ScalarActions"""

    def setUp(self):
        x = np.arange(100, dtype=float)
        x[31] = np.nan
        x[41] = np.nan
        x[59] = np.nan
        y = x**2
        self.data = {
            "r_y": x,
            "i_y": y,
            "z_y": y.reshape((10, 10)),
        }
        mask = np.zeros(100, dtype=bool)
        mask[1:50] = 1
        self.mask = {
            "r": mask,
            "i": mask,
            "z": mask.reshape((10, 10)),
        }

    def _testScalarActionAlmostEqual(self, cls, truth, maskedTruth):
        action = cls(vectorKey="{band}_y")
        schema = [col for col, colType in action.getInputSchema()]
        self.assertEqual(schema, ["{band}_y"])

        result = action(self.data, band="i")
        self.assertAlmostEqual(result, truth)

        result = action(self.data, band="z")
        self.assertAlmostEqual(result, truth)

        result = action(self.data, mask=self.mask["i"], band="i")
        self.assertAlmostEqual(result, maskedTruth)

        result = action(self.data, mask=self.mask["z"], band="z")
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
        r = np.arange(float(self.size)) + 1
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
        selector = RangeSelector(vectorKey="r_vector", minimum=0, maximum=self.size)
        prefix = "a_"
        stats = CalcBinnedStatsAction(name_prefix=prefix, selector_range=selector, key_vector="r_vector")
        result = stats(self.data)
        mask = selector(self.data)
        values = self.data["r_vector"][mask]
        median = np.median(values)
        truth = {
            stats.name_mask: mask,
            stats.name_median: median,
            stats.name_sigmaMad: 1.482602218505602 * np.median(np.abs(values - median)),
            stats.name_count: np.sum(mask),
            stats.name_select_maximum: np.max(values),
            stats.name_select_median: median,
            stats.name_select_minimum: np.min(values),
            "range_maximum": self.size,
            "range_minimum": 0,
        }
        self.assertEqual(list(result.keys()), list(truth.keys()))

        np.testing.assert_array_almost_equal(result[stats.name_sigmaMad], truth[stats.name_sigmaMad])
        del truth[stats.name_sigmaMad]

        np.testing.assert_array_equal(result[stats.name_mask], truth[stats.name_mask])
        del truth[stats.name_mask]

        for key, value in truth.items():
            self.assertEqual(result[key], value, key)

    def testCalcMomentSize(self):
        xx = self.data["r_ixx"]
        yy = self.data["r_iyy"]
        xy = self.data["r_ixy"]

        # Test determinant with defaults
        action = CalcMomentSize()
        result = action(self.data, band="r")
        schema = [col for col, colType in action.getInputSchema()]
        self.assertEqual(sorted(schema), ["{band}_ixx", "{band}_ixy", "{band}_iyy"])
        truth = 0.25 * (xx * yy - xy**2)
        np.testing.assert_array_almost_equal(result, truth)

        # Test trace with columns specified
        action = CalcMomentSize(
            colXx="{band}_iixx",
            colYy="{band}_iiyy",
            colXy="{band}_iixy",
            sizeType="trace",
        )
        result = action(self.data, band="g")
        schema = [col for col, colType in action.getInputSchema()]
        self.assertEqual(sorted(schema), ["{band}_iixx", "{band}_iiyy"])
        truth = sqrt(0.5 * (xx + yy))
        np.testing.assert_array_almost_equal(result, truth)

        CalcMomentSize(sizeType="trace", colXy=None).validate()
        with self.assertRaises(FieldValidationError):
            CalcMomentSize(sizeType="determinant", colXy=None).validate()

    # def testCalcE(self): TODO: implement

    # def testCalcEDiff(self): TODO: implement

    # def testCalcE1(self): TODO: implement

    # def testCalcE2(self): TODO: implement

    # MathActions

    def _testMath(self, ActionType, truth, num_vectors: int = 2, compare_exact: bool = False, **kwargs):
        letters = ("A", "B")
        bands = ("r", "i")
        actions = {
            f"action{letters[num]}": LoadVector(vectorKey=f"{{band{num+1}}}_vector")
            for num in range(num_vectors)
        }
        action = ActionType(**actions, **kwargs)
        kwargs_bands = {f"band{num+1}": bands[num] for num in range(num_vectors)}
        result = action(self.data, **kwargs_bands)
        self._checkSchema(action, [action.vectorKey for action in actions.values()])
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
        self._testMath(AddVector, truth, compare_exact=True)

    def testSubtract(self):
        truth = [0.0, -2.0, -6.0, -12.0, -20.0]
        self._testMath(SubtractVector, truth, compare_exact=True)

    def testMultiply(self):
        truth = [1.0, 8.0, 27.0, 64.0, 125.0]
        self._testMath(MultiplyVector, truth, compare_exact=False)

    def testDivide(self):
        truth = 1 / np.arange(1, 6)
        self._testMath(DivideVector, truth, compare_exact=False)

    def testSqrt(self):
        truth = sqrt(np.arange(1, 6))
        self._testMath(SqrtVector, truth, compare_exact=False, num_vectors=1)

    def testSquare(self):
        truth = np.arange(1, 6) ** 2
        self._testMath(SquareVector, truth, compare_exact=True, num_vectors=1)

    def testRaiseFromBase(self):
        power = np.arange(1, 6)
        for base in (-2.3, 0.6):
            truth = base**power
            self._testMath(RaiseFromBaseVector, truth, compare_exact=False, base=base, num_vectors=1)

    def testRaiseToPower(self):
        base = np.arange(1, 6)
        for power in (-2.3, 0.6):
            truth = base**power
            self._testMath(RaiseToPowerVector, truth, compare_exact=False, power=power, num_vectors=1)

    def testLog10(self):
        truth = np.log10(np.arange(1, 6))
        self._testMath(Log10Vector, truth, compare_exact=False, num_vectors=1)

    def testFractionalDifference(self):
        truth = [0.0, -0.5, -0.6666666666666666, -0.75, -0.8]
        self._testMath(FractionalDifference, truth, compare_exact=False)

    # Basic vectorActions

    # def testLoadVector(self): TODO: implement

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


class TestVectorRhoStats(unittest.TestCase):
    """Test Rho stats"""

    def setUp(self):
        # generate data just for testCalcRhoStatistics.
        np.random.seed(42)
        sizeRho = 1000
        size_src = np.random.normal(scale=1e-3, size=sizeRho)
        e1_src = np.random.normal(scale=1e-3, size=sizeRho)
        e2_src = np.random.normal(scale=1e-3, size=sizeRho)

        size_psf = np.random.normal(scale=1e-3, size=sizeRho)
        e1_psf = np.random.normal(scale=1e-3, size=sizeRho)
        e2_psf = np.random.normal(scale=1e-3, size=sizeRho)

        src_data = np.array(
            [self.getMatrixElements(size, e1, e2) for size, e1, e2 in zip(size_src, e1_src, e2_src)]
        )
        psf_data = np.array(
            [self.getMatrixElements(size, e1, e2) for size, e1, e2 in zip(size_psf, e1_psf, e2_psf)]
        )

        dataRhoStats = {
            "coord_ra": np.random.uniform(-120, 120, sizeRho),
            "coord_dec": np.random.uniform(-120, 120, sizeRho),
            "r_ixx": src_data[:, 0],
            "r_iyy": src_data[:, 1],
            "r_ixy": src_data[:, 2],
            "r_ixxPSF": psf_data[:, 0],
            "r_iyyPSF": psf_data[:, 1],
            "r_ixyPSF": psf_data[:, 2],
        }

        self.dataRhoStats = pd.DataFrame.from_dict(dataRhoStats)

    # Needed for testCalcRhoStatistics.
    @staticmethod
    def getMatrixElements(size, e1, e2):
        # puting guards just in case e1 or e2 are
        # supprior to 1, but unlikely.
        if abs(e1) >= 1:
            e1 = 0
        if abs(e2) >= 1:
            e2 = 0
        e = sqrt(e1**2 + e2**2)
        q = (1 - e) / (1 + e)
        phi = 0.5 * np.arctan2(e2, e1)
        rot = np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])
        ell = np.array([[size**2, 0], [0, (size * q) ** 2]])
        L = np.dot(rot.T, ell.dot(rot))
        return [L[0, 0], L[1, 1], L[0, 1]]

    def testCalcRhoStatistics(self):

        # just check if runs
        rho = CalcRhoStatistics()
        rho.treecorr.nbins = 21
        rho.treecorr.min_sep = 0.01
        rho.treecorr.max_sep = 100.0
        rho.treecorr.sep_units = "arcmin"
        result = rho(self.dataRhoStats, band="r")
        for rho in result:
            if rho != "rho3alt":
                self.assertEqual(np.sum(np.isfinite(result[rho].xip)), len(result[rho].xip))
                self.assertEqual(np.sum(np.isfinite(result[rho].xim)), len(result[rho].xim))
            else:
                self.assertEqual(np.sum(np.isfinite(result[rho].xi)), len(result[rho].xi))


class TestVectorSelectors(unittest.TestCase):
    def setUp(self):
        self.size = 20
        falseFlags = {
            "{band}_psfFlux_flag": [1],
            "{band}_pixelFlags_saturatedCenter": [3],
            "{band}_extendedness_flag": [5],
            "coord_flag": [7],
            "i_pixelFlags_edge": [13],
            "r_pixelFlags_edge": [15],
            "i_pixelFlags_nodata": [14],
            "r_pixelFlags_nodata": [16],
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
            "sky_object",
            "coord_flag",
            "detect_isPatchInner",
            "detect_isDeblendedSource",
        ]
        self._checkSchema(selector, keys)

        result = selector(self.data)
        truth = np.ones(self.size, dtype=bool)
        for bit in (1, 3, 5, 7, 9, 11, 13, 15, 17):
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
        selector = RangeSelector(vectorKey="r_psfFlux", minimum=np.nextafter(20, 30), maximum=50)
        self._checkSchema(selector, ["r_psfFlux"])
        result = self.data["r_psfFlux"][selector(self.data)]
        truth = [30, 40]
        np.testing.assert_array_equal(result, truth)

    def testSetSelector(self):
        n_values = 3
        values = self.data["r_psfFlux"][:n_values]
        selector = SetSelector(vectorKeys=("r_psfFlux", "i_cmodelFlux"), values=values)
        self._checkSchema(selector, ("r_psfFlux", "i_cmodelFlux"))
        result = selector(self.data)
        truth = np.zeros_like(result)
        truth[:n_values] = True
        # i_cModelFlux is just r_psfFlux reversed
        truth[-n_values:] = True
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
        keys = ["{band}_pixelFlags_edge", "{band}_pixelFlags_nodata", "sky_object"]
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


class TestKeyedDataActions(unittest.TestCase):
    def testCalcRelativeDistances(self):
        # To test CalcRelativeDistances, make a matched visit catalog with
        # objects in a box slightly larger than the annulus used in calculating
        # relative distances.
        num_visits = 15
        scatter_in_degrees = (5 * u.milliarcsecond).to(u.degree).value
        obj_id = 0
        visit_id = range(num_visits)
        all_ras, all_decs, all_objs, all_visits = [], [], [], []
        for ra in np.linspace(0, 6, 10):
            for dec in np.linspace(0, 6, 10):
                ra_degrees = (ra * u.arcmin).to(u.degree).value
                dec_degrees = (dec * u.arcmin).to(u.degree).value
                ra_meas = ra_degrees + np.random.rand(num_visits) * scatter_in_degrees
                dec_meas = dec_degrees + np.random.rand(num_visits) * scatter_in_degrees
                all_ras.append(ra_meas)
                all_decs.append(dec_meas)
                all_objs.append(np.ones(num_visits) * obj_id)
                all_visits.append(visit_id)
                obj_id += 1
        data = pd.DataFrame(
            {
                "coord_ra": np.concatenate(all_ras),
                "coord_dec": np.concatenate(all_decs),
                "obj_index": np.concatenate(all_objs),
                "visit": np.concatenate(all_visits),
            }
        )

        task = CalcRelativeDistances()
        res = task(data)

        self.assertNotEqual(res["AMx"], np.nan)
        self.assertNotEqual(res["ADx"], np.nan)
        self.assertNotEqual(res["AFx"], np.nan)


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
