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

__all__ = ("PropertyMapTool",)

from typing import Any, Dict, Iterable, MutableMapping

from healsparse.healSparseMap import HealSparseMap
from lsst.pex.config import Field

from ..actions.plot.propertyMapPlot import PropertyMapPlot
from ..interfaces import AnalysisAction, AnalysisTool


class LoadHealSparseMap(AnalysisAction):
    """Load a selection of HealSparseMaps configured for plotting."""

    mapsKey = Field[str](
        doc="The key used to access the dictionary of requested HealSparseMap objects to be loaded."
    )

    def getInputSchema(self) -> Iterable[tuple[str, Dict]]:
        return [(self.mapsKey, Dict[str, HealSparseMap])]

    def __call__(
        self, data: MutableMapping[str, Dict[str, HealSparseMap]], **kwds: Any
    ) -> Dict[str, HealSparseMap]:
        return data[self.mapsKey]


class PropertyMapTool(AnalysisTool):
    """An `AnalysisTool` for plotting property maps."""

    # Make the getOutputNames() method in the plot action config-aware.
    dynamicOutputNames: bool = True

    # Do not iterate over multiple bands in a parameterized manner.
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.data = LoadHealSparseMap()
        self.process.buildActions.data.mapsKey = "maps"
        self.produce.plot = PropertyMapPlot()
