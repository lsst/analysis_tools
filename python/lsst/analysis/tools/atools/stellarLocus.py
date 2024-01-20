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
    MagSelector,
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

        self.prep.selectors.magSelector = MagSelector()
        self.prep.selectors.magSelector.maxMag = 22.5

        self.prep.selectors.starSelector1 = StarSelector()
        self.prep.selectors.starSelector2 = StarSelector()
        self.prep.selectors.starSelector3 = StarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.returnMillimags = False

        self.process.calculateActions.approxMagDepth = ApproxFloor(vectorKey="mag")

        self.produce.plot = ColorColorFitPlot()


class WPerpPSF(StellarLocusBase):
    def setDefaults(self):
        super().setDefaults()
        colorBands = ["g", "r", "i"]
        fluxType = "psfFlux"
        self.prep.selectors.flagSelector.bands = colorBands
        self.prep.selectors.snSelector.bands = ["r"]
        self.prep.selectors.snSelector.fluxType = "{band}_" + fluxType
        self.prep.selectors.magSelector.bands = ["r"]
        self.prep.selectors.magSelector.fluxType = "{band}_" + fluxType

        self.prep.selectors.starSelector1.vectorKey = colorBands[0] + "_extendedness"
        self.prep.selectors.starSelector2.vectorKey = colorBands[1] + "_extendedness"
        self.prep.selectors.starSelector3.vectorKey = colorBands[2] + "_extendedness"

        self.process.buildActions.x.magDiff.col1 = colorBands[0] + "_" + fluxType
        self.process.buildActions.x.magDiff.col2 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col1 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col2 = colorBands[2] + "_" + fluxType
        self.process.buildActions.mag = ConvertFluxToMag(vectorKey=colorBands[1] + "_" + fluxType)

        self.process.calculateActions.wPerp_psfFlux = StellarLocusFitAction()
        stellarLocusFitDict = {
            "xMin": 0.28,
            "xMax": 1.0,
            "yMin": 0.02,
            "yMax": 0.48,
            # The "fixed" values were roughly derived from w_2023_50 HSC-RC2
            # processing (DM-42194). See commit message for further details.
            "mFixed": 0.483,
            "bFixed": -0.042,
        }
        for item in stellarLocusFitDict.items():
            self.process.calculateActions.wPerp_psfFlux.stellarLocusFitDict.update([item])

        self.produce.metric.units = {
            "wPerp_psfFlux_sigmaMAD": "mmag",
            "wPerp_psfFlux_median": "mmag",
        }

        self.produce.plot.xAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[0], colorBands[1], fluxType.replace("Flux", "")
        )
        self.produce.plot.yAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[1], colorBands[2], fluxType.replace("Flux", "")
        )
        self.produce.plot.xLims = (-0.7, 1.99)
        self.produce.plot.yLims = (-0.7, 2.49)
        self.produce.plot.magLabel = "{} {}_Mag".format(
            self.prep.selectors.magSelector.bands[0], fluxType.replace("Flux", "")
        )
        self.produce.plot.plotName = "wPerp_" + fluxType


class WPerpCModel(WPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        colorBands = ["g", "r", "i"]
        fluxType = "cModelFlux"

        self.prep.selectors.snSelector.fluxType = "{band}_" + fluxType
        self.prep.selectors.magSelector.fluxType = "{band}_" + fluxType

        self.process.buildActions.x.magDiff.col1 = colorBands[0] + "_" + fluxType
        self.process.buildActions.x.magDiff.col2 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col1 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col2 = colorBands[2] + "_" + fluxType
        self.process.buildActions.mag = ConvertFluxToMag(vectorKey=colorBands[1] + "_" + fluxType)

        self.process.calculateActions.wPerp_cModelFlux = StellarLocusFitAction()
        stellarLocusFitDict = {
            "xMin": 0.28,
            "xMax": 1.0,
            "yMin": 0.02,
            "yMax": 0.48,
            # The "fixed" values were roughly derived from w_2023_50 HSC-RC2
            # processing (DM-42194). See commit message for further details.
            "mFixed": 0.483,
            "bFixed": -0.042,
        }
        for item in stellarLocusFitDict.items():
            self.process.calculateActions.wPerp_cModelFlux.stellarLocusFitDict.update([item])

        self.produce.metric.units = {
            "wPerp_cModelFlux_sigmaMAD": "mmag",
            "wPerp_cModelFlux_median": "mmag",
        }

        self.produce.plot.xAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[0], colorBands[1], fluxType.replace("Flux", "")
        )
        self.produce.plot.yAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[1], colorBands[2], fluxType.replace("Flux", "")
        )
        self.produce.plot.magLabel = "{} {}_Mag".format(
            self.prep.selectors.magSelector.bands[0], fluxType.replace("Flux", "")
        )
        self.produce.plot.plotName = "wPerp_" + fluxType


class XPerpPSF(WPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        fluxType = "psfFlux"

        self.process.calculateActions.xPerp_psfFlux = StellarLocusFitAction()
        stellarLocusFitDict = {
            "xMin": 1.05,
            "xMax": 1.55,
            "yMin": 0.78,
            "yMax": 1.62,
            # The "fixed" values were roughly derived from w_2023_50 HSC-RC2
            # processing (DM-42194). See commit message for further details.
            "mFixed": 60.0,
            "bFixed": -75.0,
        }
        for item in stellarLocusFitDict.items():
            self.process.calculateActions.xPerp_psfFlux.stellarLocusFitDict.update([item])

        self.produce.metric.units = {
            "xPerp_psfFlux_sigmaMAD": "mmag",
            "xPerp_psfFlux_median": "mmag",
        }
        self.produce.plot.plotName = "xPerp_" + fluxType


class XPerpCModel(WPerpCModel):
    def setDefaults(self):
        super().setDefaults()
        fluxType = "cModelFlux"

        self.process.calculateActions.xPerp_cModelFlux = StellarLocusFitAction()
        stellarLocusFitDict = {
            "xMin": 1.05,
            "xMax": 1.55,
            "yMin": 0.78,
            "yMax": 1.62,
            # The "fixed" values were roughly derived from w_2023_50 HSC-RC2
            # processing (DM-42194). See commit message for further details.
            "mFixed": 60.0,
            "bFixed": -75.0,
        }
        for item in stellarLocusFitDict.items():
            self.process.calculateActions.xPerp_cModelFlux.stellarLocusFitDict.update([item])

        self.produce.metric.units = {
            "xPerp_cModelFlux_sigmaMAD": "mmag",
            "xPerp_cModelFlux_median": "mmag",
        }
        self.produce.plot.plotName = "xPerp_" + fluxType


class YPerpPSF(StellarLocusBase):
    def setDefaults(self):
        super().setDefaults()

        colorBands = ["r", "i", "z"]
        fluxType = "psfFlux"
        self.prep.selectors.flagSelector.bands = colorBands
        self.prep.selectors.snSelector.bands = ["i"]
        self.prep.selectors.snSelector.fluxType = "{band}_" + fluxType
        self.prep.selectors.magSelector.bands = ["i"]
        self.prep.selectors.magSelector.fluxType = "{band}_" + fluxType

        self.prep.selectors.starSelector1.vectorKey = colorBands[0] + "_extendedness"
        self.prep.selectors.starSelector2.vectorKey = colorBands[1] + "_extendedness"
        self.prep.selectors.starSelector3.vectorKey = colorBands[2] + "_extendedness"

        self.process.buildActions.x.magDiff.col1 = colorBands[0] + "_" + fluxType
        self.process.buildActions.x.magDiff.col2 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col1 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col2 = colorBands[2] + "_" + fluxType
        self.process.buildActions.mag = ConvertFluxToMag(vectorKey=colorBands[1] + "_" + fluxType)

        self.process.calculateActions.yPerp_psfFlux = StellarLocusFitAction()
        stellarLocusFitDict = {
            "xMin": 0.82,
            "xMax": 2.01,
            "yMin": 0.37,
            "yMax": 0.90,
            # The "fixed" values were roughly derived from w_2023_50 HSC-RC2
            # processing (DM-42194). See commit message for further details.
            "mFixed": 0.385,
            "bFixed": 0.064,
        }
        for item in stellarLocusFitDict.items():
            self.process.calculateActions.yPerp_psfFlux.stellarLocusFitDict.update([item])

        self.produce.metric.units = {
            "yPerp_psfFlux_sigmaMAD": "mmag",
            "yPerp_psfFlux_median": "mmag",
        }

        self.produce.plot.xAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[0], colorBands[1], fluxType.replace("Flux", "")
        )
        self.produce.plot.yAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[1], colorBands[2], fluxType.replace("Flux", "")
        )
        self.produce.plot.xLims = (-0.8, 3.59)
        self.produce.plot.yLims = (-0.5, 1.49)
        self.produce.plot.magLabel = "{} {}_Mag".format(
            self.prep.selectors.magSelector.bands[0], fluxType.replace("Flux", "")
        )
        self.produce.plot.plotName = "yPerp_psfFlux"


class YPerpCModel(YPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        colorBands = ["r", "i", "z"]
        fluxType = "cModelFlux"

        self.prep.selectors.snSelector.fluxType = "{band}_" + fluxType
        self.prep.selectors.magSelector.fluxType = "{band}_" + fluxType

        self.process.buildActions.x.magDiff.col1 = colorBands[0] + "_" + fluxType
        self.process.buildActions.x.magDiff.col2 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col1 = colorBands[1] + "_" + fluxType
        self.process.buildActions.y.magDiff.col2 = colorBands[2] + "_" + fluxType
        self.process.buildActions.mag = ConvertFluxToMag(vectorKey=colorBands[1] + "_" + fluxType)

        self.process.calculateActions.yPerp_cModelFlux = StellarLocusFitAction()
        stellarLocusFitDict = {
            "xMin": 0.82,
            "xMax": 2.01,
            "yMin": 0.37,
            "yMax": 0.90,
            # The "fixed" values were roughly derived from w_2023_50 HSC-RC2
            # processing (DM-42194). See commit message for further details.
            "mFixed": 0.385,
            "bFixed": 0.064,
        }
        for item in stellarLocusFitDict.items():
            self.process.calculateActions.yPerp_cModelFlux.stellarLocusFitDict.update([item])

        self.produce.metric.units = {
            "yPerp_cModelFlux_sigmaMAD": "mmag",
            "yPerp_cModelFlux_median": "mmag",
        }

        self.produce.plot.xAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[0], colorBands[1], fluxType.replace("Flux", "")
        )
        self.produce.plot.yAxisLabel = "{} - {} ({}) [mags]".format(
            colorBands[1], colorBands[2], fluxType.replace("Flux", "")
        )
        self.produce.plot.magLabel = "{} {}_Mag".format(
            self.prep.selectors.magSelector.bands[0], fluxType.replace("Flux", "")
        )
        self.produce.plot.plotName = "yPerp_" + fluxType
