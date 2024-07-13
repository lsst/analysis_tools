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

__all__ = ("IsrMetadataExposureDetectorAnalysisConfig", "IsrMetadataExposureDetectorAnalysisTask")

import pandas as pd
from lsst.pex.config import Field
from lsst.pipe.base import NoWorkFound, connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class IsrMetadataExposureDetectorAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "exposure", "detector"),
):
    metadata = connectionTypes.Input(
        doc="Task metadata from ISR",
        name="isr_metadata",
        storageClass="TaskMetadata",
        deferLoad=True,
        dimensions=["instrument", "exposure", "detector"],
    )


class IsrMetadataExposureDetectorAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=IsrMetadataExposureDetectorAnalysisConnections,
):
    subTaskName = Field[str](
        doc="The name of ISR subtask to extract metadata from. If None, the entire metadata will be used.",
        default=None,
    )


class IsrMetadataExposureDetectorAnalysisTask(AnalysisPipelineTask):
    ConfigClass = IsrMetadataExposureDetectorAnalysisConfig
    _DefaultName = "isrMetadataAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataId = butlerQC.quantum.dataId
        inputs = butlerQC.get(inputRefs)
        taskName = inputRefs.metadata.datasetType.name
        taskName = taskName[: taskName.find("_")]
        metadata = inputs["metadata"].get()
        if self.config.subTaskName:
            taskFullName = f"{taskName}:{self.config.subTaskName}"
            metadata = metadata.metadata[taskFullName].to_dict()
        else:
            taskFullName = taskName
        if not metadata:
            raise NoWorkFound(f"No metadata entries for {taskFullName}.")
        df = pd.DataFrame(metadata)

        plotInfo = self.parsePlotInfo({"data": inputs["metadata"]}, dataId)
        inputs["data"] = df
        inputs.pop("metadata")
        outputs = self.run(plotInfo=plotInfo, **inputs)
        butlerQC.put(outputs, outputRefs)
