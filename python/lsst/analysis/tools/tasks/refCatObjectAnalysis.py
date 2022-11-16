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
from ..contexts import CoaddContext
from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class RefCatObjectAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={"outputName": "objectTable_tract_gaia_dr2_20200414_match"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="objectTable_tract_gaia_dr2_20200414_match",
        storageClass="DataFrame",
        deferLoad=True,
        dimensions=("skymap", "tract"),
    )


class RefCatObjectAnalysisConfig(AnalysisBaseConfig, pipelineConnections=RefCatObjectAnalysisConnections):
    def setDefaults(self):
        super().setDefaults()

        # set plots to run
        self.plots.astromDiffRA = TargetRefCatDeltaRAScatterPlot()
        self.plots.astromDiffRA.parameterizedBand = True
        self.plots.astromDiffRA.applyContext(CoaddContext)

        self.plots.astromDiffDec = TargetRefCatDeltaDecScatterPlot()
        self.plots.astromDiffDec.parameterizedBand = True
        self.plots.astromDiffDec.applyContext(CoaddContext)

        # set metrics to run - none so far


class RefCatObjectAnalysisTask(AnalysisPipelineTask):
    """Make plots and metrics using a table of objects matched to reference
    catalog sources.
    """

    ConfigClass = RefCatObjectAnalysisConfig
    _DefaultName = "refCatObjectAnalysisTask"
