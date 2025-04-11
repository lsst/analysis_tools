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
    "MetadataAnalysisConfig",
    "DatasetMetadataAnalysisTask",
    "TaskMetadataAnalysisTask",
)

from lsst.pex.config import Field, ListField
from lsst.pipe.base import UpstreamFailureNoWorkFound, connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask
from astropy.utils.iers import conf
conf.auto_max_age = None

class MetadataAnalysisConnections(
    AnalysisBaseConnections,
    dimensions={},
    defaultTemplates={"inputName": "", "outputName": "", "datasetType": "", "storageClass": ""},
):
    def __init__(self, *, config=None):
        """Customize the connections for a specific task's or dataset's
        metadata.

        Parameters
        ----------
        config : `MetadataAnalysisConfig`
            A config for analyzing task or dataset metadata with this
            connection.
        """
        # The following must come before super().__init__ so the input
        # dimensions are propagated into the output metric bundle by
        # the base class __init__
        self.dimensions = frozenset(config.inputDimensions)
        super().__init__(config=config)

        self.data = connectionTypes.Input(
            doc="Input dataset to extract metadata from.",
            name=config.connections.inputName,
            storageClass=config.connections.storageClass,
            deferLoad=True,
            dimensions=self.dimensions,
        )


class MetadataAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=MetadataAnalysisConnections,
):
    inputDimensions = ListField(
        # Sort to ensure default order is consistent between runs
        default=sorted(MetadataAnalysisConnections.dimensions),
        dtype=str,
        doc="The dimensions of the input dataset.",
    )
    raiseNoWorkFoundOnEmptyMetadata = Field(
        dtype=bool,
        default=False,
        doc="Raise a NoWorkFound error if none of the configured metrics are in the task metadata.",
    )
    raiseNoWorkFoundOnIncompleteMetadata = Field(
        dtype=bool,
        default=False,
        doc="Raise NoWorkFound if any of the configured metrics are not in the task metadata.",
    )

    def validate(self):
        super().validate()
        if self.raiseNoWorkFoundOnEmptyMetadata and self.raiseNoWorkFoundOnIncompleteMetadata:
            raise ValueError("At most one 'raiseNoWorkFound' option may be True at a time.")


class DatasetMetadataAnalysisConfig(MetadataAnalysisConfig):

    def setDefaults(self):
        super().setDefaults()
        self.raiseNoWorkFoundOnIncompleteMetadata = True


class DatasetMetadataAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DatasetMetadataAnalysisConfig
    _DefaultName = "datasetMetadataAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataId = butlerQC.quantum.dataId
        inputs = butlerQC.get(inputRefs)
        plotInfo = self.parsePlotInfo(inputs, dataId)
        data = inputs["data"].get()

        # Collect all the metrics that are configured to run.
        metadata = {}
        for name in self.config.atools.fieldNames:
            atool = getattr(self.config.atools, name)
            if hasattr(atool, "metrics"):
                for metric in atool.metrics:
                    # Check if the metric uses base key prefixes.
                    if atool.metricsPrefixedWithBaseKeys.get(metric):
                        # Find all keys prefixed with the base key and collect
                        # them in a dictionary.
                        metadata[metric] = {k: v for k, v in data.metadata.items() if k.startswith(metric)}
                    else:
                        # Retrieve the metric directly if it's not prefixed.
                        metadata[metric] = data.metadata.get(metric)
                    # Check if the retrieved metadata is empty.
                    if not metadata[metric] and self.config.raiseNoWorkFoundOnIncompleteMetadata:
                        raise UpstreamFailureNoWorkFound(
                            f"Metadata entry '{metric}' is empty for {inputRefs.data.datasetType.name}, "
                            f"or it is not one of {data.metadata.getOrderedNames()}."
                        )
                if self.config.raiseNoWorkFoundOnEmptyMetadata and not any(v for v in metadata.values()):
                    raise UpstreamFailureNoWorkFound(
                        "All configured metadata entries are missing from "
                        f"{inputRefs.data.datasetType.name}."
                    )

        outputs = self.run(data={"metadata_metrics": metadata}, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)


class TaskMetadataAnalysisTask(AnalysisPipelineTask):
    ConfigClass = MetadataAnalysisConfig
    _DefaultName = "taskMetadataAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataId = butlerQC.quantum.dataId
        inputs = butlerQC.get(inputRefs)
        plotInfo = self.parsePlotInfo(inputs, dataId)
        metadata = inputs["data"].get().to_dict()
        taskName = inputRefs.data.datasetType.name
        taskName = taskName[: taskName.find("_")]
        if not metadata:
            raise UpstreamFailureNoWorkFound(f"No metadata entries for {taskName}.")
        if self.config.raiseNoWorkFoundOnEmptyMetadata or self.config.raiseNoWorkFoundOnIncompleteMetadata:
            self.validateMetrics(metadata, taskName)
        outputs = self.run(data=metadata, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)

    def validateMetrics(self, metadata, taskName):
        """Raise NoWorkFound if there are insufficent metrics in the task
        metadata.

        Parameters
        ----------
        metadata : `dict`
            The task metadata converted to a dict.
        taskName : `str`
            The name of the task to extract metadata from

        Raises
        ------
        NoWorkFound
            If none of the metrics are in the metadata.
        """
        for fieldName in self.config.atools.fieldNames:
            atool = getattr(self.config.atools, fieldName)
            subTaskNames = getattr(atool, "subTaskNames", None) or {}
            for key in atool.metrics.keys():
                if key in metadata[taskName].keys():
                    if self.config.raiseNoWorkFoundOnEmptyMetadata:
                        return
                elif subtaskName := subTaskNames.get(key):
                    if f"{taskName}:{subtaskName}" not in metadata:
                        raise UpstreamFailureNoWorkFound(
                            f"Subtask {subtaskName!r} was not found in {taskName} metadata"
                        )
                    if key in metadata[f"{taskName}:{subtaskName}"]:
                        if self.config.raiseNoWorkFoundOnEmptyMetadata:
                            return
                else:
                    if self.config.raiseNoWorkFoundOnIncompleteMetadata:
                        raise UpstreamFailureNoWorkFound(
                            f"Metric {key!r} was not found in {taskName} metadata"
                        )
        if self.config.raiseNoWorkFoundOnEmptyMetadata:
            raise UpstreamFailureNoWorkFound(
                f"None of the specified metrics were found in the {taskName} metadata"
            )
