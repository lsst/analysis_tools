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

__all__ = ("DiffimDetectorVisitMetricsAnalysisConfig", "DiffimDetectorVisitMetricsAnalysisTask")

import pandas as pd
from lsst.pipe.base import NoWorkFound, connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiffimDetectorVisitMetricsAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
):
    metadataSubtract = connectionTypes.Input(
        doc="Task metadata from image differencing",
        name="subtractImages_metadata",
        storageClass="TaskMetadata",
        dimensions=("visit", "band", "detector"),
    )
    metadataDetect = connectionTypes.Input(
        doc="Task metadata from DIA detection and measurement",
        name="detectAndMeasure_metadata",
        storageClass="TaskMetadata",
        dimensions=("visit", "band", "detector"),
    )


class DiffimDetectorVisitMetricsAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiffimDetectorVisitMetricsAnalysisConnections
):
    pass


class DiffimDetectorVisitMetricsAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiffimDetectorVisitMetricsAnalysisConfig
    _DefaultName = "DiffimDetectorVisitMetricsAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        detectTaskName = inputRefs.metadataDetect.datasetType.name
        detectTaskName = detectTaskName[: detectTaskName.find("_")]
        metadata = inputs["metadataDetect"].metadata[detectTaskName].to_dict()
        if not metadata:
            raise NoWorkFound("No metadata entries for detectAndMeasure.")
        inputs.pop("metadataDetect")
        subtractTaskName = inputRefs.metadataSubtract.datasetType.name
        subtractTaskName = subtractTaskName[: subtractTaskName.find("_")]
        metadata |= inputs["metadataSubtract"].metadata[subtractTaskName].to_dict()
        inputs.pop("metadataSubtract")
        df = pd.DataFrame(metadata)

        inputs["data"] = df
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)
