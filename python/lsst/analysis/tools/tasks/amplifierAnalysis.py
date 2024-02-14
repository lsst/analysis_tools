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

__all__ = ("AmplifierAnalysisConfig", "AmplifierAnalysisTask")

from lsst.pipe.base import connectionTypes as ct
from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class AmplifierAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument",),
    defaultTemplates={
        "inputDataType": "verifyBiasResults",
        "outputName": "biasPercentiles",
    },
):
    data = ct.Input(
        doc="Exposure and detector based amplifier bias distributions.",
        name="{inputDataType}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument",),
    )


class AmplifierAnalysisConfig(AnalysisBaseConfig, pipelineConnections=AmplifierAnalysisConnections):
    pass


class AmplifierAnalysisTask(AnalysisPipelineTask):
    """Make plots and metrics using tables of bias, flat, and dark
    distributions.
    """

    ConfigClass = AmplifierAnalysisConfig
    _DefaultName = "amplifierAnalysisTask"
