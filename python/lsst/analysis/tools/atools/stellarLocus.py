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

__all__ = (
    "StellarLocusBase",
    "WPerpPSF",
    "WPerpCModel",
    "XPerpPSF",
    "XPerpCModel",
    "YPerpPSF",
    "YPerpCModel",
)

from ..actions.keyedData.stellarLocusFit import StellarLocusFitAction
from ..actions.plot.colorColorFitPlot import ColorColorFitPlot
from ..actions.scalar import ApproxFloor
from ..actions.vector import (
    CoaddPlotFlagSelector,
    ConvertFluxToMag,
    ExtinctionCorrectedMagDiff,
    SnSelector,
    StarSelector,
)
from ..interfaces import AnalysisTool


class StellarLocusBase(AnalysisTool):
    # Use this as the Base Class for now StellarLocusBaseMetric
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.threshold = 300

        self.prep.selectors.starSelector = StarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.returnMillimags = False


class WPerpPSF(StellarLocusBase):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]
        self.prep.selectors.snSelector.bands = ["r"]
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"

        self.prep.selectors.starSelector.vectorKey = "r_extendedness"

        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"

        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"

        self.process.buildActions.mag = ConvertFluxToMag(vectorKey="r_psfFlux")

        self.process.calculateActions.wPerp_psfFlux = StellarLocusFitAction()
        self.process.calculateActions.wPerp_psfFlux.stellarLocusFitDict = {
            "xMin": 0.28,
            "xMax": 1.0,
            "yMin": 0.02,
            "yMax": 0.48,
            "mHW": 0.52,
            "bHW": -0.08,
        }

        self.process.calculateActions.approxMagDepth = ApproxFloor(vectorKey="mag")
        self.produce.metric.units = {
            "wPerp_psfFlux_sigmaMAD": "mmag",
            "wPerp_psfFlux_median": "mmag",
        }

        self.produce.plot = ColorColorFitPlot()
        self.produce.plot.xAxisLabel = "g - r (PSF) [mags]"
        self.produce.plot.yAxisLabel = "r - i (PSF) [mags]"
        self.produce.plot.magLabel = "PSF Mag"
        self.produce.plot.plotName = "wPerp_psfFlux"


class WPerpCModel(WPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.snSelector.fluxType = "{band}_cModelFlux"

        self.process.buildActions.x.magDiff.col1 = "g_cModelFlux"
        self.process.buildActions.x.magDiff.col2 = "r_cModelFlux"

        self.process.buildActions.y.magDiff.col1 = "r_cModelFlux"
        self.process.buildActions.y.magDiff.col2 = "i_cModelFlux"

        self.process.buildActions.mag = ConvertFluxToMag(vectorKey="r_cModelFlux")

        self.process.calculateActions.wPerp_cmodelFlux = StellarLocusFitAction()
        self.process.calculateActions.wPerp_cmodelFlux.stellarLocusFitDict = {
            "xMin": 0.28,
            "xMax": 1.0,
            "yMin": 0.02,
            "yMax": 0.48,
            "mHW": 0.52,
            "bHW": -0.08,
        }

        self.process.calculateActions.approxMagDepth = ApproxFloor(vectorKey="mag")

        self.produce.metric.units = {
            "wPerp_cmodelFlux_sigmaMAD": "mmag",
            "wPerp_cmodelFlux_median": "mmag",
        }

        self.produce.plot = ColorColorFitPlot()
        self.produce.plot.xAxisLabel = "g - r (CModel) [mags]"
        self.produce.plot.yAxisLabel = "r - i (CModel) [mags]"
        self.produce.plot.magLabel = "CModel Mag"
        self.produce.plot.plotName = "wPerp_cmodelFlux"


class XPerpPSF(StellarLocusBase):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]
        self.prep.selectors.snSelector.bands = ["r"]
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"

        self.prep.selectors.starSelector.vectorKey = "r_extendedness"

        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"

        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"

        self.process.buildActions.mag = ConvertFluxToMag(vectorKey="r_psfFlux")

        self.process.calculateActions.xPerp_psfFlux = StellarLocusFitAction()
        self.process.calculateActions.xPerp_psfFlux.stellarLocusFitDict = {
            "xMin": 1.05,
            "xMax": 1.55,
            "yMin": 0.78,
            "yMax": 1.62,
            "mHW": 13.35,
            "bHW": -15.54,
        }

        self.process.calculateActions.approxMagDepth = ApproxFloor(vectorKey="mag")

        self.produce.metric.units = {
            "xPerp_psfFlux_sigmaMAD": "mmag",
            "xPerp_psfFlux_median": "mmag",
        }
        self.produce.plot = ColorColorFitPlot()
        self.produce.plot.xAxisLabel = "g - r (PSF) [mags]"
        self.produce.plot.yAxisLabel = "r - i (PSF) [mags]"
        self.produce.plot.magLabel = "PSF Mag"
        self.produce.plot.plotName = "xPerp_psfFlux"


class XPerpCModel(XPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.snSelector.fluxType = "{band}_cModelFlux"

        self.process.buildActions.x.magDiff.col1 = "g_cModelFlux"
        self.process.buildActions.x.magDiff.col2 = "r_cModelFlux"

        self.process.buildActions.y.magDiff.col1 = "r_cModelFlux"
        self.process.buildActions.y.magDiff.col2 = "i_cModelFlux"

        self.process.calculateActions.xPerp_cmodelFlux = StellarLocusFitAction()
        self.process.calculateActions.xPerp_cmodelFlux.stellarLocusFitDict = {
            "xMin": 1.05,
            "xMax": 1.55,
            "yMin": 0.78,
            "yMax": 1.62,
            "mHW": 13.35,
            "bHW": -15.54,
        }

        self.process.buildActions.mag = ConvertFluxToMag(vectorKey="r_cModelFlux")

        self.produce.metric.units = {  # type: ignore
            "xPerp_cmodelFlux_sigmaMAD": "mmag",
            "xPerp_cmodelFlux_median": "mmag",
        }

        self.produce.plot = ColorColorFitPlot()
        self.produce.plot.xAxisLabel = "g - r (CModel) [mags]"
        self.produce.plot.yAxisLabel = "r - i (CModel) [mags]"
        self.produce.plot.magLabel = "CModel Mag"
        self.produce.plot.plotName = "xPerp_cmodelFlux"


class YPerpPSF(StellarLocusBase):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector.bands = ["r", "i", "z"]
        self.prep.selectors.snSelector.bands = ["i"]
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"

        self.prep.selectors.starSelector.vectorKey = "i_extendedness"

        self.process.buildActions.x.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "i_psfFlux"

        self.process.buildActions.y.magDiff.col1 = "i_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "z_psfFlux"

        self.process.buildActions.mag = ConvertFluxToMag(vectorKey="i_psfFlux")

        self.process.calculateActions.yPerp_psfFlux = StellarLocusFitAction()
        self.process.calculateActions.yPerp_psfFlux.stellarLocusFitDict = {
            "xMin": 0.82,
            "xMax": 2.01,
            "yMin": 0.37,
            "yMax": 0.90,
            "mHW": 0.40,
            "bHW": 0.03,
        }

        self.process.calculateActions.approxMagDepth = ApproxFloor(vectorKey="mag")

        self.produce.metric.units = {
            "yPerp_psfFlux_sigmaMAD": "mmag",
            "yPerp_psfFlux_median": "mmag",
        }

        self.produce.plot = ColorColorFitPlot()
        self.produce.plot.xAxisLabel = "r - i (PSF) [mags]"
        self.produce.plot.yAxisLabel = "i - z (PSF) [mags]"
        self.produce.plot.magLabel = "PSF Mag"
        self.produce.plot.plotName = "yPerp_psfFlux"


class YPerpCModel(YPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.snSelector.fluxType = "{band}_cModelFlux"

        self.process.buildActions.x.magDiff.col1 = "r_cModelFlux"
        self.process.buildActions.x.magDiff.col2 = "i_cModelFlux"

        self.process.buildActions.y.magDiff.col1 = "i_cModelFlux"
        self.process.buildActions.y.magDiff.col2 = "z_cModelFlux"

        self.process.buildActions.mag = ConvertFluxToMag(vectorKey="i_cModelFlux")

        self.process.calculateActions.yPerp_cmodelFlux = StellarLocusFitAction()
        self.process.calculateActions.yPerp_cmodelFlux.stellarLocusFitDict = {
            "xMin": 0.82,
            "xMax": 2.01,
            "yMin": 0.37,
            "yMax": 0.90,
            "mHW": 0.40,
            "bHW": 0.03,
        }

        self.produce.metric.units = {  # type: ignore
            "yPerp_cmodelFlux_sigmaMAD": "mmag",
            "yPerp_cmodelFlux_median": "mmag",
        }

        self.produce.plot = ColorColorFitPlot()
        self.produce.plot.xAxisLabel = "r - i (CModel) [mags]"
        self.produce.plot.yAxisLabel = "i - z (CModel) [mags]"
        self.produce.plot.magLabel = "CModel Mag"
        self.produce.plot.plotName = "yPerp_cmodelFlux"
