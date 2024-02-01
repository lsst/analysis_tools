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
    "VerifyPtcAnalysisConfig",
    "VerifyPtcAnalysisTask",
    "VerifyCalibAnalysisConfig",
    "VerifyCalibAnalysisTask",
)

from lsst.pipe.base import connectionTypes as cT

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class VerifyCalibAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument",),
    defaultTemplates={"inputName": "verifyBiasResults"},
):
    data = cT.Input(
        doc="Table containing bias verification data to load from the butler",
        name="verifyBiasResults",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
    )

    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )


class VerifyCalibAnalysisConfig(AnalysisBaseConfig, pipelineConnections=VerifyCalibAnalysisConnections):
    pass


class VerifyCalibAnalysisTask(AnalysisPipelineTask):
    ConfigClass = VerifyCalibAnalysisConfig
    _DefaultName = "verifyCalibAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class.
        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)
        data = self.loadData(inputs["data"])
        camera = inputs["camera"]

        outputs = self.run(
            data=data,
            plotInfo=plotInfo,
            camera=camera,
        )
        butlerQC.put(outputs, outputRefs)


# Photon Transfer Curve: PTC
class VerifyPtcAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument",),
    defaultTemplates={"inputName": "verifyPtcResults"},
):
    data = cT.Input(
        doc="Table containing PTC verification data to load from the butler",
        name="verifyPtcResults",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
    )

    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )


class VerifyPtcAnalysisConfig(AnalysisBaseConfig, pipelineConnections=VerifyPtcAnalysisConnections):
    def setDefaults(self):
        super().setDefaults()


class VerifyPtcAnalysisTask(AnalysisPipelineTask):
    ConfigClass = VerifyPtcAnalysisConfig
    _DefaultName = "verifyPtcAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class.
        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)
        data = self.loadData(inputs["data"])
        camera = inputs["camera"]

        outputs = self.run(
            data=data,
            plotInfo=plotInfo,
            camera=camera,
        )
        butlerQC.put(outputs, outputRefs)
