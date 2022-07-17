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

from ..analysisMetrics.analysisMetrics import ShapeSizeFractionalMetric
from ..analysisPlots.analysisPlots import Ap12_PSF_skyPlot, ShapeSizeFractionalDiffScatter
from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class ObjectTableTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={"inputName": "objectTable_tract"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="objectTable_tract",
        storageClass="DataFrame",
        # deferLoad=True,
        dimensions=("skymap", "tract"),
    )


class ObjectTableTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=ObjectTableTractAnalysisConnections
):
    def setDefaults(self):
        super().setDefaults()
        # set plots to run
        self.plots.shapeSizeFractionalDiffScatter = ShapeSizeFractionalDiffScatter()
        self.plots.Ap12_PSF_skyPlot = Ap12_PSF_skyPlot()

        # set metrics to run
        self.metrics.shapeSizeFractionalMetric = ShapeSizeFractionalMetric()


class ObjectTableTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = ObjectTableTractAnalysisConfig
    _DefaultName = "objectTableTractAnalysisTask"
