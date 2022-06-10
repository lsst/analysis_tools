from __future__ import annotations

__all__ = ("ShapeSizeFractionalDiffScatter",)

from ..analysisParts.shapeSizeFractional import (
    ShapeSizeFractionalPrep,
    ShapeSizeFractionalProcess,
    ShapeSizeFractionalPostProcessPlot
)

from ..interfaces import AnalysisPlot


class ShapeSizeFractionalDiffScatter(AnalysisPlot):
    def setDefaults(self):
        super().setDefaults()
        self.prep = ShapeSizeFractionalPrep()
        self.process = ShapeSizeFractionalProcess()
        self.post_process = ShapeSizeFractionalPostProcessPlot()

        self.post_process.highThreshold = self.process.highSNRSelector.threshold  # type: ignore
        self.post_process.lowThreshold = self.process.lowSNRSelector.threshold  # type: ignore
