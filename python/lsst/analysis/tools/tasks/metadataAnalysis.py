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

__all__ = ("MetadataAnalysisConfig", "TaskMetadataAnalysisTask", "DictTypeDatasetMetadataAnalysisTask")

from lsst.pex.config import Field, ListField
from lsst.pipe.base import NoWorkFound, connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class MetadataAnalysisConnections(
    AnalysisBaseConnections,
    dimensions={"instrument", "exposure", "detector"},
    defaultTemplates={"taskName": "", "datasetType": "", "storageClass": ""},
):
    def __init__(self, *, config=None):
        """Customize the connections for a specific task's or dataset type's
        metadata.

        Parameters
        ----------
        config : `MetadataMetricConfig`
            A config for `MetadataMetricTask` or one of its subclasses.
        """
        super().__init__(config=config)

        if set(config.metadataDimensions) != self.dimensions:
            self.dimensions.clear()
            self.dimensions.update(config.metadataDimensions)

        if config.isTaskMetadata:
            if not config.connections.taskName:
                raise ValueError("Task name must be set to use task metadata.")
            if config.connections.datasetType or config.connections.storageClass:
                raise ValueError("Dataset type and storage class are not used with task metadata.")
            name = f"{config.connections.taskName}_metadata"
            doc = "Task metadata to load."
            storageClass = "TaskMetadata"
        else:
            if not config.connections.datasetType or not config.connections.storageClass:
                raise ValueError("Dataset type and storage class must be set to use dataset metadata.")
            if config.connections.taskName:
                raise ValueError("Task name is not used with dataset metadata.")
            name = f"{config.connections.datasetType}"
            doc = "Input dataset type to load."
            storageClass = f"{config.connections.storageClass}"

        self.data = connectionTypes.Input(
            doc=doc,
            name=name,
            storageClass=storageClass,
            deferLoad=True,
            dimensions=frozenset(config.metadataDimensions),
        )


class MetadataAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=MetadataAnalysisConnections,
):
    isTaskMetadata = Field[bool](
        doc="Whether to use task metadata instead of metadata from a dataset type.",
        default=True,
    )

    metadataDimensions = ListField(
        # Sort to ensure default order is consistent between runs
        default=sorted(MetadataAnalysisConnections.dimensions),
        dtype=str,
        doc="Override for the dimensions of the input data.",
    )


class TaskMetadataAnalysisTask(AnalysisPipelineTask):
    ConfigClass = MetadataAnalysisConfig
    _DefaultName = "metadataAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataId = butlerQC.quantum.dataId
        inputs = butlerQC.get(inputRefs)
        plotInfo = self.parsePlotInfo(inputs, dataId)
        metadata = inputs["data"].get().to_dict()
        if not metadata:
            taskName = inputRefs.data.datasetType.name
            taskName = taskName[: taskName.find("_")]
            raise NoWorkFound(f"No metadata entries for {taskName}.")
        outputs = self.run(data=metadata, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)


class DictTypeDatasetMetadataAnalysisTask(AnalysisPipelineTask):
    ConfigClass = MetadataAnalysisConfig
    _DefaultName = "dictTypeMetadataAnalysis"

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
                    metadata[metric] = data.metadata.get_dict(metric)
                    if not metadata[metric]:
                        raise NoWorkFound(
                            f"Metadata entry '{metric}' is empty for {inputRefs.data.datasetType.name}."
                        )

        outputs = self.run(data={"metadata_metrics": metadata}, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)
