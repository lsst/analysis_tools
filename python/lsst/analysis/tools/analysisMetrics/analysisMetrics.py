from __future__ import annotations
from lsst.analysis.tools.vectorActions.calcShapeSize import CalcShapeSize

from lsst.analysis.tools.vectorActions.vectorActions import (
    DownselectVector,
    ExtinctionCorrectedMagDiff,
    FractionalDifference,
    MagColumnNanoJansky,
)

from ..vectorActions.selectors import CoaddPlotFlagSelector, SnSelector, StellarSelector, VectorSelector

# from ..analysisParts.wPerp import WPerpProcess, WPerpPostProcessMetric
from ..interfaces import AnalysisMetric, KeyedData
from ..plotActions.scatterplotWithTwoHists import ScatterPlotStatsAction
from ..keyedDataActions.stellarLocusFit import StellarLocusFitAction
from typing import Any


class ShapeSizeFractionalMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()  # type: ignore
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)  # type: ignore

        self.process.buildActions.mags = MagColumnNanoJansky(columnKey="{band}_psfFlux")  # type: ignore
        self.process.buildActions.fracDiff = FractionalDifference(  # type: ignore
            actionA=CalcShapeSize(),
            actionB=CalcShapeSize(colXx="{band}_ixxPSF", colYy="{band}_iyyPSF", colXy="{band}_ixyPSF"),
        )
        # pre-compute a stellar selector mask so it can be used in the filter
        # actions while only being computed once, alternatively the stellar
        # selector could be calculated and applied twice in the filter stage
        self.process.buildActions.starSelector = StellarSelector()  # type: ignore

        self.process.filterActions.xStars = DownselectVector(  # type: ignore
            vectorKey="mags", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.yStars = DownselectVector(  # type: ignore
            vectorKey="fracDiff", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.stars = ScatterPlotStatsAction(  # type: ignore
            vectorKey="fracDiff",
        )

        self.post_process.units = {  # type: ignore
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        }


class WPerpPSFMetric(AnalysisMetric):
    # Use this as the Base Class for now StellarLocusBaseMetric
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
        self.process.calculateActions.wPerp.stellarLocusFitDict = {"xMin": 0.1, "xMax": 0.2,
                                                                   "yMin": 0.1, "yMax": 0.2,
                                                                   "mHW": 0.5, "bHW": 0.0}
        self.post_process.units = {  # type: ignore
            "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }

    def __call__(self, data: KeyedData, **kwargs) -> Any:
        kwargs.pop("bands")
        return super().__call__()


class XPerpMetric(AnalysisMetric):
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

        self.post_process.units = {  # type: ignore
            "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class WPerpPSFMetric(WPerpPSFMetric):
    pass


class WPerpCModelMetric(WPerpPSFMetric):
    pass


class XPerpPSFMetric(WPerpPSFMetric):
    pass


class XPerpCModelMetric(WPerpPSFMetric):
    pass


class YPerpPSFMetric(WPerpPSFMetric):
    pass


class YPerpCModelMetric(WPerpPSFMetric):
    pass
