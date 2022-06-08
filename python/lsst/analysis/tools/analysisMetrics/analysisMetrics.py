from __future__ import annotations

from ..analysisParts.shapeSizeFractional import (
    ShapeSizeFractionalPrep,
    ShapeSizeFractionalProcess,
    ShapeSizeFractionalPostProcessMetric,
)
from ..interfaces import AnalysisMetric


class ShapeSizeFractionalMetric(AnalysisMetric):
    def setDefaults(self):
        # configure the prep step
        self.prep = ShapeSizeFractionalPrep()
        self.process = ShapeSizeFractionalProcess()
        self.post_process = ShapeSizeFractionalPostProcessMetric()
