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
    CoaddPlotFlagSelector,
    ExtinctionCorrectedMagDiff,
    MagColumnNanoJansky,
    SnSelector,
    StarSelector,
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
