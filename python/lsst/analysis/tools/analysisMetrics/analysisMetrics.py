from __future__ import annotations

from ..analysisParts.shapeSizeFractional import (
    ShapeSizeFractionalPrep,
    ShapeSizeFractionalProcessMetric,
    ShapeSizeFractionalPostProcessMetric,
)
from ..interfaces import AnalysisMetric


class ShapeSizeFractionalMetric(AnalysisMetric):
    def setDefaults(self):
        # configure the prep step
        self.prep = ShapeSizeFractionalPrep()
        self.process = ShapeSizeFractionalProcessMetric()
        self.post_process = ShapeSizeFractionalPostProcessMetric()
