from __future__ import annotations

from ..actions.keyedData.stellarLocusFit import StellarLocusFitAction
from ..actions.scalar import ApproxFloor
from ..actions.vector import (
    CoaddPlotFlagSelector,
    ExtinctionCorrectedMagDiff,
    MagColumnNanoJansky,
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

        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="r_psfFlux")

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


class WPerpCModel(WPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.snSelector.fluxType = "{band}_cModelFlux"

        self.process.buildActions.x.magDiff.col1 = "g_cModelFlux"
        self.process.buildActions.x.magDiff.col2 = "r_cModelFlux"

        self.process.buildActions.y.magDiff.col1 = "r_cModelFlux"
        self.process.buildActions.y.magDiff.col2 = "i_cModelFlux"

        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="r_cModelFlux")

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

        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="r_psfFlux")

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

        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="r_cModelFlux")


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

        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="i_psfFlux")

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


class YPerpCModel(YPerpPSF):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.snSelector.fluxType = "{band}_cModelFlux"

        self.process.buildActions.x.magDiff.col1 = "r_cModelFlux"
        self.process.buildActions.x.magDiff.col2 = "i_cModelFlux"

        self.process.buildActions.y.magDiff.col1 = "i_cModelFlux"
        self.process.buildActions.y.magDiff.col2 = "z_cModelFlux"

        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="i_cModelFlux")

        self.process.calculateActions.yPerp_cmodelFlux = StellarLocusFitAction()
        self.process.calculateActions.yPerp_cmodelFlux.stellarLocusFitDict = {
            "xMin": 0.82,
            "xMax": 2.01,
            "yMin": 0.37,
            "yMax": 0.90,
            "mHW": 0.40,
            "bHW": 0.03,
        }
