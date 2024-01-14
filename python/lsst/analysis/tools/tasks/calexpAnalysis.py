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

__all__ = (
    "CalexpAnalysisConfig",
    "CalexpAnalysisTask",
)

from collections.abc import Iterable
from ..interfaces._interfaces import KeyedData
from lsst.pipe.base import InputQuantizedConnection, OutputQuantizedConnection, QuantumContext
from lsst.pipe.base import connectionTypes as ct
from lsst.pex.config import ListField

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CalexpAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"inputName": "calexp", "outputName": "pixelMaskMetrics"},
):
    data = ct.Input(
        doc="Calibrated exposure to load from the butler",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("visit", "band", "detector"),
        deferLoad=False,
    )


class CalexpAnalysisConfig(AnalysisBaseConfig, pipelineConnections=CalexpAnalysisConnections):
    images = ListField[str](doc="Images to extract from the exposure.", default=["image", "mask", "variance"])


class CalexpAnalysisTask(AnalysisPipelineTask):
    ConfigClass = CalexpAnalysisConfig
    _DefaultName = "calexpAnalysis"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        # Docstring inherited.

        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId

        imageDict = {}
        for image in self.config.images:
            imageDict[image] = getattr(inputs["data"], image).getArray()

        outputs = self.run(data={"imageArrays": imageDict})
