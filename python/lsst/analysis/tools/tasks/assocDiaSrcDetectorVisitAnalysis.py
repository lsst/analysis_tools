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

__all__ = ("AssocDiaSrcDetectorVisitAnalysisConfig", "AssocDiaSrcDetectorVisitAnalysisTask")

from lsst.pipe.base import connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask
from astropy.utils.iers import conf
conf.auto_max_age = None


class AssocDiaSrcDetectorVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": ""},
):
    data = connectionTypes.Input(
        doc="CcdVisit-based DiaSource table to load from the butler",
        name="{fakesType}{coaddName}Diff_assocDiaSrc",
        storageClass="DataFrame",
        deferLoad=True,
        dimensions=("visit", "band", "detector"),
    )


class AssocDiaSrcDetectorVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssocDiaSrcDetectorVisitAnalysisConnections
):
    pass


class AssocDiaSrcDetectorVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AssocDiaSrcDetectorVisitAnalysisConfig
    _DefaultName = "AssocDiaSrcDetectorVisitAnalysis"
