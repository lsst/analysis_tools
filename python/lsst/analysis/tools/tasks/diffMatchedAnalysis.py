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

from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiffMatchedAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={"inputName": "diff_matched_truth_summary_objectTable_tract"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="matched_truth_summary_objectTable_tract",
        storageClass="DataFrame",
        deferLoad=True,
        dimensions=("skymap", "tract"),
    )


class DiffMatchedAnalysisConfig(AnalysisBaseConfig, pipelineConnections=DiffMatchedAnalysisConnections):
    def setDefaults(self):
        super().setDefaults()


class DiffMatchedAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiffMatchedAnalysisConfig
    _DefaultName = "DiffMatchedAnalysisTask"
