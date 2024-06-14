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

__all__ = ("IsrDetectorExposureMetricsAnalysisConfig", "IsrAmpOffsetMetricsAnalysisTask")

import pandas as pd
from lsst.pipe.base import NoWorkFound, connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class IsrDetectorExposureMetricsAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "exposure", "detector"),
):
    metadata = connectionTypes.Input(
        doc="Task metadata from ISR",
        name="isr_metadata",
        storageClass="TaskMetadata",
        dimensions=["instrument", "exposure", "detector"],
    )


class IsrDetectorExposureMetricsAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=IsrDetectorExposureMetricsAnalysisConnections
):
    pass


class IsrAmpOffsetMetricsAnalysisTask(AnalysisPipelineTask):
    ConfigClass = IsrDetectorExposureMetricsAnalysisConfig
    _DefaultName = "isrAmpOffsetMetricsAnalysis"
    _SubTaskName = "ampOffset"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        taskName = inputRefs.metadata.datasetType.name
        taskName = taskName[: taskName.find("_")]
        metadata = inputs["metadata"].metadata[f"{taskName}:{self._SubTaskName}"].to_dict()
        if not metadata:
            raise NoWorkFound(f"No metadata entries for {taskName}:{self._SubTaskName}.")
        inputs.pop("metadata")
        df = pd.DataFrame(metadata)

        inputs["data"] = df
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)
