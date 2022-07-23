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

__all__ = [
    "StellarLocusBaseMetric",
    "WPerpPSFMetric",
]

from ..actions.keyedData.stellarLocusFit import StellarLocusFitAction
from ..actions.scalar.scalarActions import ApproxFloor
from ..actions.vector import (
    AndSelector,
    BandSelector,
    CalcShapeSize,
    CoaddPlotFlagSelector,
    ExtinctionCorrectedMagDiff,
    MagColumnNanoJansky,
    PerGroupStatistic,
    SkyObjectSelector,
    Sn,
    SnSelector,
    StarSelector,
    ThresholdSelector,
    VectorSelector,
)
from ..interfaces import AnalysisMetric


class StellarLocusBaseMetric(AnalysisMetric):
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


class GRIStellarPSFMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
    parameterizedBand: bool = False
    fluxType: str = "psfFlux"

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_" + self.fluxType
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "r_extendedness"

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = f"g_{self.fluxType}"
        self.process.buildActions.x.magDiff.col2 = f"r_{self.fluxType}"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = f"r_{self.fluxType}"
        self.process.buildActions.y.magDiff.col2 = f"i_{self.fluxType}"
        self.process.buildActions.y.magDiff.returnMillimags = False
        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="r_psfFlux")

        self.process.calculateActions.approxMagDepth = ApproxFloor(vectorKey="mag")
        self.process.calculateActions.wPerp_psfFlux = StellarLocusFitAction()
        self.process.calculateActions.wPerp_psfFlux.stellarLocusFitDict = {
            "xMin": 0.28,
            "xMax": 1.0,
            "yMin": 0.02,
            "yMax": 0.48,
            "mHW": 0.52,
            "bHW": -0.08,
        }

        # self.process.calculateActions.xPerp = StellarLocusFitAction()
        # self.process.calculateActions.xPerp.stellarLocusFitDict = {}

        self.produce.units = {  # type: ignore
            "wPerp_psfFlux_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


#            "xPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp


class GRIStellarCModelMetric(GRIStellarPSFMetric):
    fluxType: str = "CModel"

    def setDefaults(self):
        super().setDefaults()

        self.produce.newNames = {
            "wPerp_sigmaMAD": "wCmodelPerp_sigmaMAD",  # TODO need to return mmag from wPerp
            "xPerp_sigmaMAD": "xCmodelPerp_sigmaMAD",  # TODO need to return mmag from wPerp
        }


class RIZStellarPSFMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["r", "i", "z"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["i"]

        self.prep.selectors.starSelector = StarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "i_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = "i_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "z_psfFlux"
        self.process.buildActions.y.magDiff.returnMillimags = False

        self.process.calculateActions.yPerp = StellarLocusFitAction()
        self.process.calculateActions.yerp.stellarLocusFitDict = {}

        self.produce.units = {  # type: ignore
            "yPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class WPerpPSFMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"
        self.process.buildActions.y.magDiff.returnMillimags = False

        self.process.calculateActions.wPerp = StellarLocusFitAction()
        self.process.calculateActions.wPerp.stellarLocusFitDict = {
            "xMin": 0.28,
            "xMax": 1.0,
            "yMin": 0.02,
            "yMax": 0.48,
            "mHW": 0.52,
            "bHW": -0.08,
        }
        self.produce.units = {  # type: ignore
            "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class XPerpPSFMetric(AnalysisMetric):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"
        self.process.buildActions.y.magDiff.returnMillimags = False

        self.process.calculateActions.xPerp = StellarLocusFitAction()

        self.produce.units = {  # type: ignore
            "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class SkyFluxStatisticMetric(AnalysisMetric):
    parameterizedBand: bool = True
    fluxType: str = "ap09Flux"

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.skyObjectSelector = SkyObjectSelector()
        self.prep.selectors.skyObjectSelector.bands = ["{band}"]

        self.process.calculateActions.medianSky = MedianAction(vectorKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.meanSky = MeanAction(vectorKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.stdevSky = StdevAction(vectorKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.sigmaMADSky = SigmaMadAction(vectorKey=f"{{band}}_{self.fluxType}")

        self.produce.units = {  # type: ignore
            "medianSky": "nJy",
            "meanSky": "nJy",
            "stdevSky": "nJy",
            "sigmaMADSky": "nJy",
        }


class StellarPhotometricRepeatabilityMetric(AnalysisMetric):
    parameterizedBand: bool = True
    fluxType: str = "psfFlux"

    def setDefaults(self):
        super().setDefaults()

        # Apply per-source selection criteria
        self.prep.selectors.bandSelector = BandSelector()

        # Compute per-group quantities
        self.process.buildActions.perGroupSn = PerGroupStatistic()
        self.process.buildActions.perGroupSn.buildAction = Sn(fluxType=f"{self.fluxType}")
        self.process.buildActions.perGroupSn.func = "median"
        self.process.buildActions.perGroupExtendedness = PerGroupStatistic()
        self.process.buildActions.perGroupExtendedness.buildAction.vectorKey = "extendedness"
        self.process.buildActions.perGroupExtendedness.func = "median"
        self.process.buildActions.perGroupCount = PerGroupStatistic()
        self.process.buildActions.perGroupCount.buildAction.vectorKey = f"{self.fluxType}"
        self.process.buildActions.perGroupStdev = PerGroupStatistic()
        self.process.buildActions.perGroupStdev.buildAction = MagColumnNanoJansky(vectorKey=f"{self.fluxType}")
        self.process.buildActions.perGroupStdev.func = "std"

        # Filter on per-group quantities
        self.process.filterActions.perGroupStdev = DownselectVector(vectorKey="perGroupStdev")
        self.process.filterActions.perGroupStdev.selector = AndSelector()
        self.process.filterActions.perGroupStdev.selector.selectors.count = ThresholdSelector(
            columnKey="perGroupCount", op="ge", threshold=3,
        )
        self.process.filterActions.perGroupStdev.selector.selectors.sn = ThresholdSelector(
            columnKey="perGroupSn", op="ge", threshold=100,
        )
        self.process.filterActions.perGroupStdev.selector.selectors.extendedness = ThresholdSelector(
            columnKey="perGroupExtendedness", op="le", threshold=0.5,
        )

        self.process.calculateActions.photRepeatStdev = MedianAction(colKey="perGroupStdev")

        self.produce.units = {  # type: ignore
            "photRepeatStdev": "mag",
        }
        self.produce.newNames = {
            "photRepeatStdev": "{band}_stellarPhotRepeatStdev"
        }
