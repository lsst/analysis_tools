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
    "CalexpSummaryAnalysisConfig",
    "CalexpSummaryAnalysisTask",
)

import dataclasses

import lsst.pex.config as pexConfig
from lsst.pipe.base import (
    InputQuantizedConnection,
    OutputQuantizedConnection,
    QuantumContext,
    UpstreamFailureNoWorkFound,
)
from lsst.pipe.base import connectionTypes as cT

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CalexpSummaryAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"inputName": "calexp.summaryStats", "outputName": "calexpSummary"},
):
    data = cT.Input(
        doc="Calibrated exposure summary statistics to load from the butler",
        name="calexp.summaryStats",
        storageClass="ExposureSummaryStats",
        dimensions=("visit", "band", "detector"),
        deferLoad=False,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        parent = self.data.name.split(".")[0]
        if config.input_image_type == "future":
            name = f"{parent}.summary_stats"
            storage_class = "ObservationSummaryStats"
        else:
            name = f"{parent}.summaryStats"
            storage_class = "ExposureSummaryStats"
        self.data = dataclasses.replace(self.data, name=name, storageClass=storage_class)


class CalexpSummaryAnalysisConfig(AnalysisBaseConfig, pipelineConnections=CalexpSummaryAnalysisConnections):
    input_image_type = pexConfig.ChoiceField[str](
        "Which image type to read the summary statistics from.",
        allowed={
            "legacy": "Read ``summaryStats`` from an `lsst.afw.image.Exposure`.",
            "future": "Read ``summary_stats`` from an `lsst.images.VisitImage`.",
        },
        optional=False,
        default="legacy",
    )


class CalexpSummaryAnalysisTask(AnalysisPipelineTask):
    ConfigClass = CalexpSummaryAnalysisConfig
    _DefaultName = "calexpSummary"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        # Docstring inherited.

        inputs = butlerQC.get(inputRefs)

        summary = inputs["data"]
        if summary is None:
            raise UpstreamFailureNoWorkFound("No summary stats attached to the input image.")
        outputs = self.run(data=summary.__dict__)
        butlerQC.put(outputs, outputRefs)
