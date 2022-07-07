from __future__ import annotations
# from msilib.schema import Media

from ..vectorActions.vectorActions import (
    DownselectVector,
    ExtinctionCorrectedMagDiff,
    FractionalDifference,
    MagColumnNanoJansky,
)
from ..vectorActions.calcShapeSize import CalcShapeSize

from ..vectorActions.selectors import (
    CoaddPlotFlagSelector, SnSelector,
    StellarSelector, VectorSelector,
    FlagSelector, SkyObjectSelector,
)

# from ..analysisParts.wPerp import WPerpProcess, WPerpPostProcessMetric
from ..interfaces import AnalysisMetric
from ..plotActions.scatterplotWithTwoHists import ScatterPlotStatsAction
from ..keyedDataActions.stellarLocusFit import StellarLocusFitAction
from lsst.pex.config import Field
from ..scalarActions.scalarActions import (
    ApproxFloor, MedianAction, MeanAction, SigmaMadAction,
    StdevAction,
)


class ShapeSizeFractionalMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)

        self.process.buildActions.mags = MagColumnNanoJansky(columnKey="{band}_psfFlux")
        self.process.buildActions.fracDiff = FractionalDifference(
            actionA=CalcShapeSize(),
            actionB=CalcShapeSize(colXx="{band}_ixxPSF", colYy="{band}_iyyPSF", colXy="{band}_ixyPSF"),
        )
        # pre-compute a stellar selector mask so it can be used in the filter
        # actions while only being computed once, alternatively the stellar
        # selector could be calculated and applied twice in the filter stage
        self.process.buildActions.starSelector = StellarSelector()

        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="fracDiff", selector=VectorSelector(vectorKey="starSelector")
        )
        # downselect the psfFlux as well
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="{band}_psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="{band}_psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.stars = ScatterPlotStatsAction(
            vectorKey="yStars",
        )
        # use the downselected psfFlux
        self.process.calculateActions.stars.highSNSelector.fluxType = 'psfFlux'
        self.process.calculateActions.stars.lowSNSelector.fluxType = 'psfFlux'
        self.process.calculateActions.stars.fluxType = 'psfFlux'

        self.post_process.units = {  # type: ignore
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        }


class StellarLocusBaseMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
    multiband: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        # self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        # self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        # self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StellarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        # self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        # self.process.buildActions.x.magDiff.col2 = "r_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        # self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        # self.process.buildActions.y.magDiff.col2 = "i_psfFlux"
        self.process.buildActions.y.magDiff.returnMillimags = False

        # self.process.calculateActions.wPerp = StellarLocusFitAction()
        #self.process.calculateActions.wPerp.stellarLocusFitDict = {"xMin": 0.28, "xMax": 1.0,
        #                                                           "yMin": 0.02, "yMax": 0.48,
        #                                                           "mHW": 0.52, "bHW": -0.08}
        #self.post_process.units = {  # type: ignore
        #    "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        #}


class GRIStellarPSFMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
    multiband: bool = False
    fluxType: str = "psfFlux"

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_" + self.fluxType
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StellarSelector()
        self.prep.selectors.starSelector.columnKey = 'r_extendedness'

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = f"g_{self.fluxType}"
        self.process.buildActions.x.magDiff.col2 = f"r_{self.fluxType}"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = f"r_{self.fluxType}"
        self.process.buildActions.y.magDiff.col2 = f"i_{self.fluxType}"
        self.process.buildActions.y.magDiff.returnMillimags = False
        self.process.buildActions.mag = MagColumnNanoJansky(columnKey="r_psfFlux")

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

        self.post_process.units = {  # type: ignore
            "wPerp_psfFlux_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }
#            "xPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp


class GRIStellarCModelMetric(GRIStellarPSFMetric):
    fluxType: str = "CModel"

    def setDefaults(self):
        super().setDefaults()

        self.post_process.newNames = {
            "wPerp_sigmaMAD": "wCmodelPerp_sigmaMAD",  # TODO need to return mmag from wPerp
            "xPerp_sigmaMAD": "xCmodelPerp_sigmaMAD",  # TODO need to return mmag from wPerp
        }


class RIZStellarPSFMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
    multiband: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["r", "i", "z"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["i"]

        self.prep.selectors.starSelector = StellarSelector()

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

        self.post_process.units = {  # type: ignore
            "yPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class WPerpPSFMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
    multiband: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StellarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"
        self.process.buildActions.y.magDiff.returnMillimags = False

        self.process.calculateActions.wPerp = StellarLocusFitAction()
        self.process.calculateActions.wPerp.stellarLocusFitDict = {"xMin": 0.28, "xMax": 1.0,
                                                                   "yMin": 0.02, "yMax": 0.48,
                                                                   "mHW": 0.52, "bHW": -0.08}
        self.post_process.units = {  # type: ignore
            "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class XPerpPSFMetric(AnalysisMetric):
    multiband: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StellarSelector()

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"
        self.process.buildActions.y.magDiff.returnMillimags = False

        self.process.calculateActions.xPerp = StellarLocusFitAction()

        self.post_process.units = {  # type: ignore
            "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class SkyFluxStatisticMetric(AnalysisMetric):
    multiband: bool = True
    fluxType: str = "ap09Flux"

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.skyObjectSelector = SkyObjectSelector()
        self.prep.selectors.skyObjectSelector.bands = ["{band}"]

        self.process.calculateActions.medianSky = MedianAction(colKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.meanSky = MeanAction(colKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.stdevSky = StdevAction(colKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.sigmaMADSky = SigmaMadAction(colKey=f"{{band}}_{self.fluxType}")

        self.post_process.units = {  # type: ignore
            "medianSky": "nJy",
            "meanSky": "nJy",
            "stdevSky": "nJy",
            "sigmaMADSky": "nJy",
        }
