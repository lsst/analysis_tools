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
    "DiaFakesDetectorVisitAnalysisConfig",
    "DiaFakesDetectorVisitAnalysisTask",
    "AssocDiaFakesDetectorVisitAnalysisConfig",
    "AssocDiaFakesDetectorVisitAnalysisTask",
)

from lsst.pipe.base.connectionTypes import Input

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiaFakesDetectorVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "visit", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": "fakes_"},
):
    data = Input(
        doc="CcdVisit-based Matched fake to load from the butler",
        name="{fakesType}{coaddName}Diff_matchDiaSrc",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument", "visit", "detector"),
    )


class DiaFakesDetectorVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiaFakesDetectorVisitAnalysisConnections
):
    pass


class DiaFakesDetectorVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiaFakesDetectorVisitAnalysisConfig
    _DefaultName = "diaFakesDetectorVisitAnalysis"


class AssocDiaFakesDetectorVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "visit", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": "fakes_"},
):
    data = Input(
        doc="CcdVisit-based Matched fake to load from the butler",
        name="{fakesType}{coaddName}Diff_matchAssocDiaSrc",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument", "visit", "detector"),
    )


class AssocDiaFakesDetectorVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssocDiaFakesDetectorVisitAnalysisConnections
):
    pass


class AssocDiaFakesDetectorVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AssocDiaFakesDetectorVisitAnalysisConfig
    _DefaultName = "assocDiaFakesDetectorVisitAnalysis"
