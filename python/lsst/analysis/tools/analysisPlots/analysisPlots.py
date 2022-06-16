from __future__ import annotations

__all__ = ("ShapeSizeFractionalDiffScatter",)

from ..analysisParts.shapeSizeFractional import (
    ShapeSizeFractionalProcess,
)

from ..vectorActions.selectors import CoaddPlotFlagSelector, SnSelector, StellarSelector
from ..plotActions.scatterplotWithTwoHists import ScatterPlotWithTwoHists
from ..plotActions.colorColorFitPlot import ColorColorFitPlot
from ..vectorActions.vectorActions import ExtinctionCorrectedMagDiff
from ..keyedDataActions.stellarLocusFit import StellarLocusFitAction

from ..interfaces import AnalysisPlot


class ShapeSizeFractionalDiffScatter(AnalysisPlot):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()  # type: ignore
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)  # type: ignore

        self.process = ShapeSizeFractionalProcess()
        self.post_process = ScatterPlotWithTwoHists()

        self.post_process.plotTypes = ["stars"]
        self.post_process.xAxisLabel = "PSF Magnitude (mag)"
        self.post_process.yAxisLabel = "Fractional size residuals (S/S_PSF - 1)"
        self.post_process.magLabel = "PSF Magnitude (mag)"

        self.post_process.highThreshold = self.process.highSNRSelector.threshold  # type: ignore
        self.post_process.lowThreshold = self.process.lowSNRSelector.threshold  # type: ignore

class WPerpPSFPlot(AnalysisPlot):
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

        self.process = ColorColorFitPlot()

