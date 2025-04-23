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

__all__ = ("ColumnMagnitudeScatterPlot",)

from lsst.pex.config import Field

from ..actions.vector import DownselectVector, LoadVector, Log10Vector
from ..actions.vector.selectors import CoaddPlotFlagSelector, VectorSelector, VisitPlotFlagSelector
from .genericBuild import ExtendednessTool
from .genericProduce import MagnitudeScatterPlot


class ColumnMagnitudeScatterPlot(ExtendednessTool, MagnitudeScatterPlot):
    """A MagnitudeScatterPlot with a single column value on the y axis."""

    key_y = Field[str](default=None, doc="Key of column to plot on the y axis")
    log10_y = Field[bool](default=False, doc="Whether to plot log10 of the values on the y axis")

    def coaddContext(self) -> None:
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = []

    def visitContext(self) -> None:
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()

    def finalize(self):
        # A lazy check for whether finalize has already been called
        classes = self.get_classes()
        if hasattr(self.process.filterActions, self.get_name_attr_values(classes[0])):
            return
        super().finalize()
        if not self.key_y:
            raise ValueError("Must specify key to plot on y axis")

        action = LoadVector(vectorKey=self.key_y)
        if self.log10_y:
            action = Log10Vector(action=action)

        attr_column = f"get_{self.key_y}"
        setattr(self.process.buildActions, attr_column, action)

        classes = self.get_classes()
        for object_class in classes:
            setattr(
                self.process.filterActions,
                self.get_name_attr_values(object_class),
                DownselectVector(
                    vectorKey=attr_column,
                    selector=VectorSelector(vectorKey=self.get_name_attr_selector(object_class)),
                ),
            )

        if not self.produce.plot.yAxisLabel:
            self.produce.plot.yAxisLabel = f"log10({self.key_y})" if self.log10_y else self.key_y
