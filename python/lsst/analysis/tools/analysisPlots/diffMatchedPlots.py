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
from __future__ import annotations

__all__ = ("MatchedRefCoaddPlot", "MatchedRefCoaddCModelFluxPlot")

from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.vector.vectorActions import DownselectVector, VectorSelector
from ..analysisParts.diffMatched import MatchedRefCoaddDiffMagTool, MatchedRefCoaddDiffPositionTool
from ..interfaces import AnalysisPlot


class MatchedRefCoaddPlot(AnalysisPlot):
    def setDefaults(self):
        super().setDefaults()
        self.produce = ScatterPlotWithTwoHists()

        self.produce.plotTypes = ["galaxies", "stars"]
        self.produce.xAxisLabel = "Reference Magnitude (mag)"


class MatchedRefCoaddCModelPlot(MatchedRefCoaddPlot):
    def setDefaults(self):
        super().setDefaults()
        self.produce.magLabel = "cModel mag"

        # downselect the cModelFlux as well
        for prefix, plural in (("star", "Stars"), ("galaxy", "Galaxies")):
            for suffix in ("", "Err"):
                setattr(
                    self.process.filterActions,
                    f"{prefix}_cModelFlux{suffix}",
                    DownselectVector(
                        vectorKey=f"{{band}}_cModelFlux{suffix}",
                        selector=VectorSelector(vectorKey=f"{prefix}Selector"),
                    ),
                )

            statAction = ScatterPlotStatsAction(vectorKey=f"y{plural.capitalize()}")
            fluxType = f"{prefix}_cModelFlux"
            statAction.highSNSelector.fluxType = fluxType
            statAction.highSNSelector.threshold = 200
            statAction.lowSNSelector.fluxType = fluxType
            statAction.lowSNSelector.threshold = 10
            statAction.fluxType = fluxType
            setattr(self.process.calculateActions, plural, statAction)


class MatchedRefCoaddCModelFluxPlot(MatchedRefCoaddCModelPlot, MatchedRefCoaddDiffMagTool):
    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.produce.yAxisLabel = "cModel - Reference mmag"

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.produce.yAxisLabel = "chi = (cModel - Reference mag)/error"

    def setDefaults(self):
        super().setDefaults()


class MatchedRefCoaddPositionPlot(MatchedRefCoaddCModelPlot, MatchedRefCoaddDiffPositionTool):
    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.produce.yAxisLabel = f"{self.variable} position (pix)"

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.produce.yAxisLabel = f"chi = (slot - Reference {self.variable} position)/error"

    def setDefaults(self):
        super().setDefaults()
