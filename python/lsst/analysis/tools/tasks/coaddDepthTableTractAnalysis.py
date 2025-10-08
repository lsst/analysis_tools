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
    "CoaddDepthTableTractAnalysisConnections",
    "CoaddDepthTableTractAnalysisConfig",
    "CoaddDepthTableTractAnalysisTask",
)

import lsst.pex.config as pexConfig
from lsst.pipe.base import connectionTypes as cT

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CoaddDepthTableTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"coaddName": ""},
):
    data = cT.Input(
        doc="Table with coadd depth statistics based on n_image values.",
        name="{coaddName}_coadd_depth_table",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
        deferLoad=True,
    )


class CoaddDepthTableTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=CoaddDepthTableTractAnalysisConnections
):
    load_skymap = pexConfig.Field[bool](doc="Whether to load the skymap.", default=True)


class CoaddDepthTableTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = CoaddDepthTableTractAnalysisConfig
    _DefaultName = "coaddDepthTableTractAnalysis"
