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
    "CalexpSummaryConfig",
    "CalexpSummaryTask",
)

from lsst.pipe.base import connectionTypes as cT

from lsst.pipe.base import InputQuantizedConnection, OutputQuantizedConnection, QuantumContext
from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CalexpSummaryConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"inputName": "calexp", "outputName": "calexpSummaryMetrics"},
):
    data = cT.Input(
        doc="Calibrated exposure to load from the butler",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("visit", "band", "detector"),
        deferLoad=False,
    )


class CalexpSummaryConfig(AnalysisBaseConfig, pipelineConnections=CalexpSummaryConnections):
    pass


class CalexpSummaryTask(AnalysisPipelineTask):
    ConfigClass = CalexpSummaryConfig
    _DefaultName = "calexpSummary"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        # Docstring inherited.

        inputs = butlerQC.get(inputRefs)

        summary = inputs["data"].getInfo().getSummaryStats().__dict__

        outputs = self.run(data=summary)
        butlerQC.put(outputs, outputRefs)
