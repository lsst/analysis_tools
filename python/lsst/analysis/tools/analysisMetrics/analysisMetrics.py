from __future__ import annotations

from ..vectorActions.vectorActions import (
    DownselectVector,
    ExtinctionCorrectedMagDiff,
    FractionalDifference,
    MagColumnNanoJansky,
)
from ..vectorActions.calcShapeSize import CalcShapeSize

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
        # downselect the psfFlux as well
        self.process.filterActions.psfFlux = DownselectVector(  # type: ignore
            vectorKey="{band}_psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(  # type: ignore
            vectorKey="{band}_psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.stars = ScatterPlotStatsAction(  # type: ignore
            vectorKey="yStars",
        )
        # use the downselected psfFlux
        self.process.calculateActions.stars.highSNSelector.fluxType = 'psfFlux'  # type: ignore
        self.process.calculateActions.stars.lowSNSelector.fluxType = 'psfFlux'  # type: ignore
        self.process.calculateActions.stars.fluxType = 'psfFlux'  # type: ignore

        self.post_process.units = {  # type: ignore
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        }


class XPerpMetric(AnalysisMetric):
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

        self.post_process.units = {  # type: ignore
            "wPerp_sigmaMAD": "mag",  # TODO need to return mmag from wPerp
        }


class WPerpPSFMetric(XPerpMetric):
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
