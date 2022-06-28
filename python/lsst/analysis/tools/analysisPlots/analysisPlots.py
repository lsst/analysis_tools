from __future__ import annotations

__all__ = ("ShapeSizeFractionalDiffScatter",)

from ..vectorActions.selectors import CoaddPlotFlagSelector, SnSelector, StellarSelector
from ..plotActions.scatterplotWithTwoHists import ScatterPlotWithTwoHists, ScatterPlotStatsAction
from ..plotActions.colorColorFitPlot import ColorColorFitPlot
from ..plotActions.skyPlot import SkyPlot
from ..vectorActions.vectorActions import (
    ExtinctionCorrectedMagDiff,
    MagColumnNanoJansky,
    DownselectVector,
    VectorSelector,
    FractionalDifference,
)
from ..vectorActions.calcShapeSize import CalcShapeSize

from ..scalarActions.scalarActions import ApproxFloor
from ..keyedDataActions.stellarLocusFit import StellarLocusFitAction

from ..interfaces import AnalysisPlot


class ShapeSizeFractionalDiffScatter(AnalysisPlot):
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

        self.post_process = ScatterPlotWithTwoHists()

        self.post_process.plotTypes = ["stars"]
        self.post_process.xAxisLabel = "PSF Magnitude (mag)"
        self.post_process.yAxisLabel = "Fractional size residuals (S/S_PSF - 1)"
        self.post_process.magLabel = "PSF Magnitude (mag)"


class WPerpPSFPlot(AnalysisPlot):
    # WPerp does not support running in multiband mode
    multiband: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 300
        self.prep.selectors.snSelector.bands = ["r"]

        self.prep.selectors.starSelector = StellarSelector()
        self.prep.selectors.starSelector.columnKey = "r_extendedness"

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"
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

        self.post_process = ColorColorFitPlot()
        self.post_process.plotName = "wPerp_psfFlux"

class Ap12_PSF_skyPlot(AnalysisPlot):

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["{band}"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 300

        self.prep.selectors.starSelector = StellarSelector()
        self.prep.selectors.starSelector.columnKey = "{band}_extendedness"

        self.process.buildActions.z = ExtinctionCorrectedMagDiff()
        self.process.buildActions.z.magDiff.col1 = "{band}_ap12Flux"
        self.process.buildActions.z.magDiff.col2 = "{band}_psfFlux"

        self.post_process = SkyPlot()
        self.post_process.plotName = "ap12-psf_{band}"
