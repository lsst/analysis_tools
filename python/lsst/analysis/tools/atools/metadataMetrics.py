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
    "DatasetMetadataMetricTool",
    "TaskMetadataMetricTool",
)


from lsst.pex.config import DictField, Field

from ..actions.keyedData import KeyedDataKeyAccessAction
from ..interfaces import AnalysisTool


class MetadataMetricTool(AnalysisTool):
    """Base class for tools designed to extract values from metadata"""

    parameterizedBand = Field[bool](
        doc="Does this MetadataMetricTool support band as a name parameter?", default=False
    )

    metrics = DictField[str, str](doc="The metrics to extract from the metadata and their respective units.")

    newNames = DictField[str, str](
        doc="New names to allocate to the extracted metrics. Keys are the current "
        "names, values are the new names.",
        default=None,
        optional=True,
    )

    @staticmethod
    def makeValidAttributeName(name):
        """Make a valid attribute name using a simple replacement."""
        return name.replace(" ", "_")

    def validate(self):
        for metric in self.metrics.keys():
            if not self.makeValidAttributeName(metric).isidentifier():
                raise ValueError(
                    f"{metric=} must be a valid identifier after replacing spaces with underscores."
                )


class DatasetMetadataMetricTool(MetadataMetricTool):
    """Tool designed to extract values from metadata of data products"""

    metricsPrefixedWithBaseKeys = DictField[str, bool](
        doc="Whether metrics are prefixed with base keys. For each key (a metric name or base key), "
        "the corresponding boolean value specifies if the metric should be extracted using the base "
        "key as a prefix.",
        default={},
        optional=True,
    )

    def finalize(self):
        for metric in self.metrics.keys():
            setattr(
                self.process.filterActions,
                self.makeValidAttributeName(metric),
                KeyedDataKeyAccessAction(topLevelKey="metadata_metrics"),
            )
        self.produce.metric.units = dict(self.metrics.items())

        if self.newNames is not None:
            self.produce.metric.newNames = dict(self.newNames.items())


class TaskMetadataMetricTool(MetadataMetricTool):
    """This tool is designed to extract values from task metadata"""

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

    def finalize(self):
        for metric in self.metrics.keys():
            if self.subTaskNames is not None and metric in self.subTaskNames:
                taskFullName = f"{self.taskName}:{self.subTaskNames[metric]}"
            else:
                taskFullName = self.taskName
            setattr(
                self.process.filterActions,
                self.makeValidAttributeName(metric),
                KeyedDataKeyAccessAction(topLevelKey=taskFullName),
            )
        self.produce.metric.units = dict(self.metrics.items())

        if self.newNames is not None:
            self.produce.metric.newNames = dict(self.newNames.items())
