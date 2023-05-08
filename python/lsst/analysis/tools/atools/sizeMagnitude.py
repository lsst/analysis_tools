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

from ..actions.vector import CoaddPlotFlagSelector, VisitPlotFlagSelector
from .genericBuild import ExtendednessTool, SizeTool
from .genericProduce import MagnitudeScatterPlot


class SizeMagnitudePlot(ExtendednessTool, SizeTool, MagnitudeScatterPlot):
    def coaddContext(self) -> None:
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = []

    def visitContext(self) -> None:
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()

    def finalize(self):
        # TODO: Investigate why MagnitudeScatterPlot.finalize(self) is called
        # always, even if super().finalize() is omitted
        super().finalize()
        if not self.produce.yAxisLabel:
            size = self.sizes[self.size_y]
            self.produce.yAxisLabel = (
                f"log10({size.name_size}/{size.unit_size})"
                if size.log10_size
                else f"{size.name_size} ({size.unit_size})"
            )
