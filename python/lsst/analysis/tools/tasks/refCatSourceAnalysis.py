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

from ..analysisPlots.analysisPlots import TargetRefCatDeltaDecScatterPlot, TargetRefCatDeltaRAScatterPlot
from ..contexts import VisitContext
from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class RefCatSourceAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "skymap"),
    defaultTemplates={"outputName": "sourceTable_visit_astrometryRefCat_match"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="sourceTable_visit_astrometryRefCat_match",
        storageClass="DataFrame",
        deferLoad=True,
        dimensions=("visit",),
    )


class RefCatSourceAnalysisConfig(AnalysisBaseConfig, pipelineConnections=RefCatSourceAnalysisConnections):
    def setDefaults(self):
        super().setDefaults()

        # set plots to run
        self.plots.astromDiffRA = TargetRefCatDeltaRAScatterPlot()
        self.plots.astromDiffRA.applyContext(VisitContext)

        self.plots.astromDiffDec = TargetRefCatDeltaDecScatterPlot()
        self.plots.astromDiffDec.applyContext(VisitContext)

        # set metrics to run - none so far


class RefCatSourceAnalysisTask(AnalysisPipelineTask):
    """Make plots and metrics using a table of objects matched to reference
    catalog sources.
    """

    ConfigClass = RefCatSourceAnalysisConfig
    _DefaultName = "refCatSourceAnalysisTask"
