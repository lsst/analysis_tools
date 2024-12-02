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

__all__ = ("MetadataMetricTool", "TaskMetadataMetricTool", )


from lsst.pex.config import DictField, Field

from ..actions.keyedData import KeyedDataKeyAccessAction
from ..interfaces import AnalysisTool


class MetadataMetricTool(AnalysisTool):
    """This tool is designed to extract values from metadata of data products"""

    parameterizedBand = Field[bool](
        doc="Does this MetadataMetricTool support band as a name parameter?", default=False
    )

    metrics = DictField[str, str](
        doc="The metrics to extract from the metadata and their respective units."
    )

    newNames = DictField[str, str](
        doc="New names to allocate to the extracted metrics. Keys are the current "
        "names, values are the new names.",
        default=None,
        optional=True,
    )

    def finalize(self):
        for metric, unit in self.metrics.items():
            validMetricName = metric.replace(" ", "_")
            setattr(
                self.process.filterActions, f"{validMetricName}", KeyedDataKeyAccessAction(topLevelKey="metadata_metrics")
            )
        self.produce.metric.units = dict(self.metrics.items())

        if self.newNames is not None:
            self.produce.metric.newNames = dict(self.newNames.items())


class TaskMetadataMetricTool(AnalysisTool):
    """This tool is designed to extract values from task metadata"""

    parameterizedBand = Field[bool](
        doc="Does this MetadataMetricTool support band as a name parameter?", default=False
    )

    taskName = Field[str](
        doc="The name of the task to extract metadata from.",
        default=None,
    )

    subTaskNames = DictField[str, str](
        doc="The names of the subtasks to extract metadata from. "
        "If the metric name is identified as one of the keys, then "
        "the corresponding value is taken as the subTask metadata "
        "from which to extract metadata.",
        default=None,
        optional=True,
    )

    metrics = DictField[str, str](
        doc="The metrics to extract from the task metadata and their respective units."
    )

    newNames = DictField[str, str](
        doc="New names to allocate to the extracted metrics. Keys are the current "
        "names, values are the new names.",
        default=None,
        optional=True,
    )

    def finalize(self):
        for metric, unit in self.metrics.items():
            if self.subTaskNames is not None and metric in self.subTaskNames:
                taskFullName = f"{self.taskName}:{self.subTaskNames[metric]}"
            else:
                taskFullName = self.taskName
            setattr(
                self.process.filterActions, f"{metric}", KeyedDataKeyAccessAction(topLevelKey=taskFullName)
            )
        self.produce.metric.units = dict(self.metrics.items())

        if self.newNames is not None:
            self.produce.metric.newNames = dict(self.newNames.items())
