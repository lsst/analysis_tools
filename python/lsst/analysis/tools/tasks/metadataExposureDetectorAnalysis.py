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

__all__ = ("MetadataExposureDetectorAnalysisConfig", "MetadataExposureDetectorAnalysisTask")

from lsst.pipe.base import NoWorkFound, connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class MetadataExposureDetectorAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "exposure", "detector"),
    defaultTemplates={"inputName": "isr_metadata", "outputName": "isr_metadata_exposure_detector_analysis"},
):
    data = connectionTypes.Input(
        doc="Task metadata to load.",
        name="{inputName}",
        storageClass="TaskMetadata",
        deferLoad=True,
        dimensions=["instrument", "exposure", "detector"],
    )


class MetadataExposureDetectorAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=MetadataExposureDetectorAnalysisConnections,
):
    pass


class MetadataExposureDetectorAnalysisTask(AnalysisPipelineTask):
    ConfigClass = MetadataExposureDetectorAnalysisConfig
    _DefaultName = "metadataExposureDetectorAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataId = butlerQC.quantum.dataId
        inputs = butlerQC.get(inputRefs)
        plotInfo = self.parsePlotInfo(inputs, dataId)
        taskName = inputRefs.data.datasetType.name
        taskName = taskName[: taskName.find("_")]
        metadata = inputs["data"].get().to_dict()
        if not metadata:
            raise NoWorkFound(f"No metadata entries for {taskName}.")
        outputs = self.run(data=metadata, plotInfo=plotInfo)
        self.putByBand(butlerQC, outputs, outputRefs)
