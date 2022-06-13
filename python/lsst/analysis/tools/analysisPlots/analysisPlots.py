from __future__ import annotations

__all__ = ("ShapeSizeFractionalDiffScatter",)

from ..analysisParts.shapeSizeFractional import (
    ShapeSizeFractionalProcess,
)

from ..vectorActions.selectors import CoaddPlotFlagSelector, SnSelector
from ..plotActions.scatterplotWithTwoHists import ScatterPlotWithTwoHists

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
