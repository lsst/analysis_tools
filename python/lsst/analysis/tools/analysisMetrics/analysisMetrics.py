from __future__ import annotations
from lsst.analysis.tools.vectorActions.calcShapeSize import CalcShapeSize

from lsst.analysis.tools.vectorActions.vectorActions import (
    DownselectVector,
    FractionalDifference,
    MagColumnNanoJansky,
)


from ..vectorActions.selectors import CoaddPlotFlagSelector, SnSelector, StellarSelector, VectorSelector
from ..analysisParts.shapeSizeFractional import (
    ShapeSizeFractionalProcess,
)

# from ..analysisParts.wPerp import WPerpProcess, WPerpPostProcessMetric
# in the future this will be
# from ..analysisMetrics.base import BaseMetric
from ..interfaces import AnalysisMetric
from ..plotActions.scatterplotWithTwoHists import ScatterPlotStatsAction


class ShapeSizeFractionalMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()  # type: ignore
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)  # type: ignore

        self.process.buildActions.mags = MagColumnNanoJansky(columnKey="{band}_psfFlux")
        self.process.buildActions.fracDiff = FractionalDifference(
            actionA=CalcShapeSize(),
            actionB=CalcShapeSize(colXx="{band}_ixxPSF", colYy="{band}_iyyPSF", colXy="{band}_ixyPSF"),
        )
        # pre-compute a stellar selector mask so it can be used in the filter
        # actions while only being computed once, alternatively the stellar
        # selector could be calculated and applied twice in the filter stage
        self.process.buildActions.starSelector = StellarSelector()

        self.process.filterActions.xsStars = DownselectVector(
            vectorKey="mags", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.ysStars = DownselectVector(
            vectorKey="fracDiff", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.metrics = ScatterPlotStatsAction(
            vectorKey="fracDiff",
        )

        self.post_process.units = {
                "{band}_highSNStars_median": "pixel",
                "{band}_highSNStars_sigmaMad": "pixel",
                "{band}_highSNStars_count": "count",
                "{band}_lowSNStars_median": "pixel",
                "{band}_lowSNStars_sigmaMad": "pixel",
                "{band}_lowSNStars_count": "count",
        }


class WPerpMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()
        self.process = WPerpProcess()
