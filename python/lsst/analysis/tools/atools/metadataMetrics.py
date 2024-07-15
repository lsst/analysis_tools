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

__all__ = ("MetadataMetricTool",)

from typing import Any, Dict, Iterable

from lsst.pex.config import DictField, Field

from ..interfaces import AnalysisAction, AnalysisTool


class LoadDoubleKeyedData(AnalysisAction):
    """Load double keyed data"""

    name = Field[str](doc="The name of the key to load.")

    def getInputSchema(self) -> Iterable[tuple[str, Dict]]:
        return [(self.name, Dict[str, float])]

    def __call__(self, data, **kwds: Any) -> Dict[str, Dict[str, float]]:
        return data


class MetadataMetricTool(AnalysisTool):
    """This tool is designed to extract values from task metadata"""

    metrics = DictField[str, str](
        doc="The metrics to extract from the task metadata and their respective units."
    )

    def finalize(self):
        for name, unit in self.metrics.items():
            setattr(self.process.buildActions, f"{name}", LoadDoubleKeyedData(name=name))
        self.produce.metric.units = dict(self.metrics.items())
