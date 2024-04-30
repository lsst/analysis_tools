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
    "VerifyCalibAnalysisConfig",
    "VerifyCalibAnalysisTask",
    "VerifyCalibAnalysisConfigByFilter",
    "VerifyCalibAnalysisTaskByFilter",
    "VerifyCalibDetectorConfig",
    "VerifyCalibDetectorTask",
    "VerifyCalibDetectorConfigByFilter",
    "VerifyCalibDetectorTaskByFilter",
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


class VerifyCalibAnalysisConnectionsByFilter(
    AnalysisBaseConnections,
    dimensions=("instrument", "physical_filter"),
    defaultTemplates={"inputName": "verifyFlatResults"},
):
    data = cT.Input(
        doc="Table containing bias verification data to load from the butler",
        name="verifyBiasResults",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "physical_filter"),
        deferLoad=True,
    )

    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )


class VerifyCalibAnalysisConfigByFilter(
    AnalysisBaseConfig, pipelineConnections=VerifyCalibAnalysisConnectionsByFilter
):
    pass


class VerifyCalibAnalysisTaskByFilter(VerifyCalibAnalysisTask):
    ConfigClass = VerifyCalibAnalysisConfigByFilter
    _DefaultName = "verifyCalibAnalysisByFilter"

    pass


class VerifyCalibDetectorConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "detector"),
    defaultTemplates={"inputName": "verifyFlatResults"},
):
    data = cT.Input(
        doc="Table containing bias verification data to load from the butler",
        name="verifyBiasDetResults",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "detector"),
        deferLoad=True,
    )

    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )


class VerifyCalibDetectorConfig(AnalysisBaseConfig, pipelineConnections=VerifyCalibDetectorConnections):
    pass


class VerifyCalibDetectorTask(VerifyCalibAnalysisTask):
    ConfigClass = VerifyCalibDetectorConfig
    _DefaultName = "verifyCalibDetector"

    pass


class VerifyCalibDetectorConnectionsByFilter(
    AnalysisBaseConnections,
    dimensions=("instrument", "detector", "physical_filter"),
    defaultTemplates={"inputName": "verifyFlatDetResults"},
):
    data = cT.Input(
        doc="Table containing bias verification data to load from the butler",
        name="verifyBiasDetResults",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "detector", "physical_filter"),
        deferLoad=True,
    )

    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )


class VerifyCalibDetectorConfigByFilter(
    AnalysisBaseConfig, pipelineConnections=VerifyCalibDetectorConnectionsByFilter
):
    pass


class VerifyCalibDetectorTaskByFilter(VerifyCalibDetectorTask):
    ConfigClass = VerifyCalibDetectorConfigByFilter
    _DefaultName = "verifyCalibDetectorByFilter"

    pass
