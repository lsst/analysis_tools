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

__all__ = ("StructPlotAction",)

from typing import Iterable

from lsst.pex.config import Config
from lsst.pex.config.configurableActions import ConfigurableActionStructField

from ..interfaces import KeyedData, KeyedDataSchema, PlotAction, PlotResultType


class StructPlotAction(PlotAction):
    """A struct of named PlotActions."""

    actions = ConfigurableActionStructField[str, PlotAction](doc="Named plot actions")

    def getInputSchema(self) -> KeyedDataSchema:
        for action in self.actions:
            yield from action.getInputSchema()

    def getOutputNames(self, config: Config | None = None) -> Iterable[str]:
        names = []
        for key, action in self.actions.items():
            names_action = action.getOutputNames(config=config)
            if not names_action:
                names_action = [key]
            names.extend(names_action)
        return tuple(names)

    def __call__(self, data: KeyedData, **kwargs) -> PlotResultType:
        results = {}
        for key, action in self.actions.items():
            result = action(data=data, **kwargs)
            results[key] = result
        return results
