# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ("SkyBrightnessPrecisionHistPlot",)

from ..actions.plot.histPlot import HistPanel, HistPlot, HistStatsPanel
from ..actions.scalar import FracInRange, MedianAction, MaxAction
from ..actions.keyedData import KeyedScalars
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class SkyBrightnessPrecisionHistPlot(AnalysisTool):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.hist_all = LoadVector(vectorKey="SBRatio")

        self.process.calculateActions.frac_1pct = FracInRange(
            vectorKey="SBRatio",
            minimum=0.99,          
            maximum=1.01,
            percent=False,         
        )
        self.process.calculateActions.median = MedianAction(vectorKey="SBRatio")
        self.process.calculateActions.maximum = MaxAction(vectorKey="SBRatio")

        self.produce.plot = HistPlot()

        p = HistPanel()
        p.label = "All SBRatio"
        p.hists = dict(hist_all="All SBRatio")
        p.rangeType = "sigmaMad"
        p.lowerRange = 3.5
        p.upperRange = 3.5
        p.referenceValue = 1.0

        p.statsPanel = HistStatsPanel()
        p.statsPanel.statsLabels = ["Frac in [0.99,1.01]", "Median", "Max"]
        p.statsPanel.stat1 = ["frac_1pct"]  
        p.statsPanel.stat2 = ["median"]
        p.statsPanel.stat3 = ["maximum"]


        p.validate()
        self.produce.plot.panels["panel_all"] = p
