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


from lsst.pex.config import DictField, Field

from ..actions.keyedData import KeyedDataKeyAccessAction
from ..interfaces import AnalysisTool


class MetadataMetricTool(AnalysisTool):
    """This tool is designed to extract values from task metadata"""

    parameterizedBand = Field[bool](
        doc="Does this MetadataMetricTool support band as a name parameter?", default=False
    )

    taskName = Field[str](
        doc="The name of the task to extract metadata from.",
        default=None,
    )

    subTaskName = Field[str](
        doc="The name of the subtask to extract metadata from. "
        "If None, the entire task metadata will be used.",
        default=None,
    )

    metrics = DictField[str, str](
        doc="The metrics to extract from the task metadata and their respective units."
    )

    def finalize(self):
        if self.subTaskName:
            taskFullName = f"{self.taskName}:{self.subTaskName}"
        else:
            taskFullName = self.taskName
        for metric, unit in self.metrics.items():
            setattr(
                self.process.filterActions, f"{metric}", KeyedDataKeyAccessAction(topLevelKey=taskFullName)
            )
        self.produce.metric.units = dict(self.metrics.items())
