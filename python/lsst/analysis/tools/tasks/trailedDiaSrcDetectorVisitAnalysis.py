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

__all__ = ("TrailedDiaSrcDetectorVisitAnalysisConfig", "TrailedDiaSrcDetectorVisitAnalysisTask")

from lsst.pipe.base import connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask
from astropy.utils.iers import conf
conf.auto_max_age = None

class TrailedDiaSrcDetectorVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": ""},
):
    data = connectionTypes.Input(
        doc="Output temporarily storing long trailed diaSources.",
        name="{fakesType}{coaddName}Diff_longTrailedSrc",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("visit", "band", "detector"),
    )


class TrailedDiaSrcDetectorVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=TrailedDiaSrcDetectorVisitAnalysisConnections
):
    pass


class TrailedDiaSrcDetectorVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = TrailedDiaSrcDetectorVisitAnalysisConfig
    _DefaultName = "TrailedDiaSrcDetectorVisitAnalysis"
