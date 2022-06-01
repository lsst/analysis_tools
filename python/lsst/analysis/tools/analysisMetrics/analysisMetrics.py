from __future__ import annotations

from ..analysisParts.shapeSizeFractional import ShapeSizeFractionalPrep
from ..interfaces import AnalysisMetric
from ..tabularActions import TableOfScalars


class ShapeSizeFractionalMetric(AnalysisMetric):
    def setDefaults(self):
        # configure the prep step
        self.prep = ShapeSizeFractionalPrep()

        self.process = TableOfScalars()
