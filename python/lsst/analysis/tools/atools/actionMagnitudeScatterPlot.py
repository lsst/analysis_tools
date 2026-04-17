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

__all__ = ("ActionMagnitudeScatterPlot", "SingleColumnMagnitudeScatterPlot")

import logging

from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ..actions.vector import DownselectVector, LoadVector
from ..actions.vector.selectors import CoaddPlotFlagSelector, VectorSelector, VisitPlotFlagSelector
from ..interfaces import VectorAction
from .genericProduce import MagnitudeScatterPlot

_LOG = logging.getLogger(__name__)


class ActionMagnitudeScatterPlot(MagnitudeScatterPlot):
    """A MagnitudeScatterPlot with a single column value on the y axis."""

    action_vector = ConfigurableActionField[VectorAction](
        default=None,
        doc="The VectorAction returning data to plot",
    )
    key_y = Field[str](default=None, doc="Key of derived column to plot on the y axis")

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
        self.validate()

        attr_column = f"get_{self.key_y}"
        setattr(self.process.buildActions, attr_column, self.action_vector)

        for object_class in classes:
            setattr(
                self.process.filterActions,
                self.get_name_attr_values(object_class),
                DownselectVector(
                    vectorKey=attr_column,
                    selector=VectorSelector(vectorKey=self.get_name_attr_selector(object_class)),
                ),
            )


class SingleColumnMagnitudeScatterPlot(ActionMagnitudeScatterPlot):
    """A magnitude scatter plot loading a single column."""

    _default_vectorKey = "__SingleColumnMagnitudeScatterPlot__default"

    def setDefaults(self):
        super().setDefaults()
        # Unfortunately, this can't be left None because validate will get
        # called before finalize
        self.action_vector = LoadVector(vectorKey=self._default_vectorKey)

    def finalize(self):
        if not isinstance(self.action_vector, LoadVector):
            _LOG.warning(
                f"{self.action_vector=} should be an instance of LoadVector but is not; if it has no "
                f"vectorKey attribute, this action may fail."
            )
        if self.action_vector.vectorKey == self._default_vectorKey:
            self.action_vector.vectorKey = self.key_y
        super().finalize()
