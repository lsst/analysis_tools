from __future__ import annotations
from ..analysisParts.base import BasePrep
from ..analysisParts.shapeSizeFractional import (
    ShapeSizeFractionalPrep,
    ShapeSizeFractionalProcess,
    ShapeSizeFractionalPostProcessMetric,
)
from ..analysisParts.wPerp import WPerpProcess, WPerpPostProcessMetric
# in the future this will be
# from ..analysisMetrics.base import BaseMetric
from ..interfaces import AnalysisMetric


class ShapeSizeFractionalMetric(AnalysisMetric):
    def setDefaults(self):
        self.prep = ShapeSizeFractionalPrep()
        self.process = ShapeSizeFractionalProcess()
        self.post_process = ShapeSizeFractionalPostProcessMetric()


class WPerpMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()
        self.process = WPerpProcess()
        self.post_process = WPerpPostProcessMetric()
