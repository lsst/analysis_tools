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

__all__ = (
    "PerTractPropertyMapTool",
    "SurveyWidePropertyMapTool",
)

from lsst.pex.config import Field

from ..actions.healSparseMap.healSparseMapActions import LoadHealSparseMap
from ..actions.plot.propertyMapPlot import PerTractPropertyMapPlot, SurveyWidePropertyMapPlot
from ..interfaces import AnalysisTool


class PerTractPropertyMapTool(AnalysisTool):
    """An `AnalysisTool` for plotting per-tract property maps."""

    # Do not iterate over multiple bands in a parameterized manner.
    parameterizedBand: bool = False

    nBinsHist = Field(
        dtype=int,
        doc="Number of bins to use for the histogram.",
        default=100,
    )

    def finalize(self):
        super().finalize()
        self.process.buildActions.data.mapKey = self._name.split(".")[1]

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.data = LoadHealSparseMap(mapKey="")
        self.produce.plot = PerTractPropertyMapPlot()
        self.produce.plot.plotName = "Per-Tract HealSparse Property Map"


class SurveyWidePropertyMapTool(AnalysisTool):
    """An `AnalysisTool` for plotting survey-wide property maps."""

    # Do not iterate over multiple bands in a parameterized manner.
    parameterizedBand: bool = False

    def finalize(self):
        super().finalize()
        self.process.buildActions.data.mapKey = self._name.split(".")[1]

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.data = LoadHealSparseMap(mapKey="")
        self.produce.plot = SurveyWidePropertyMapPlot()
        self.produce.plot.plotName = "Survey-Wide HealSparse Property Map"
