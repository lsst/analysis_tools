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
    "VerifyBiasResultsAnalysisConnections",
    "VerifyBiasResultsAnalysisConfig",
    "VerifyBiasResultsAnalysisTask",
)

from lsst.pipe.base.connectionTypes import Input

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class VerifyBiasResultsAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument",),
):
    data = Input(
        doc="verifyBiasResults",
        name="verifyBiasResults",
        storageClass="ArrowAstropy",
        dimensions=("instrument","day_obs"),
        deferLoad=True,
    )


class VerifyBiasResultsAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=VerifyBiasResultsAnalysisConnections,
):
    pass


class VerifyBiasResultsAnalysisTask(AnalysisPipelineTask):
    ConfigClass = VerifyBiasResultsAnalysisConfig
    _DefaultName = "verifyBiasResultsAnalysis"
