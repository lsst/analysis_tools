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

__all__ = ("DiffimDetectorVisitAnalysisConfig", "DiffimDetectorVisitAnalysisTask")

import pandas as pd
from lsst.pipe.base import connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiffimDetectorVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
):
    metadataSubtract = ct.Input(
        doc="Task metadata from image differencing",
        name="subtractImages_metadata",
        storageClass="TaskMetadata",
        dimensions=("visit", "band", "detector"),
    )
    metadataDetect = ct.Input(
        doc="Task metadata from DIA detection and measurement",
        name="detectAndMeasure_metadata",
        storageClass="TaskMetadata",
        dimensions=("visit", "band", "detector"),
    )


class DiffimDetectorVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiffimDetectorVisitAnalysisConnections
):
    pass


class DiffimDetectorVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiffimDetectorVisitAnalysisConfig
    _DefaultName = "DiffimDetectorVisitAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        metadata = inputs["metadataDetect"].metadata["detectAndMeasure"].to_dict()
        inputs.pop("metadataDetect")
        metadata |= inputs["metadataSubtract"].metadata["subtractImages"].to_dict()
        inputs.pop("metadataSubtract")
        df = pd.DataFrame(metadata)

        inputs["data"] = df
        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)
