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

from typing import TYPE_CHECKING, Any, Mapping

__all__ = ("DiffimDetectorVisitSpatiallySampledPlotsConfig", "DiffimDetectorVisitSpatiallySampledPlotsTask")

from lsst.pipe.base import connectionTypes
from lsst.utils import inheritDoc

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

if TYPE_CHECKING:
    from lsst.daf.butler import DataCoordinate


class DiffimDetectorVisitSpatiallySampledPlotsConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": ""},
):
    data = connectionTypes.Input(
        doc="QA metrics evaluated in locations throughout the difference image.",
        name="{fakesType}{coaddName}Diff_spatiallySampledMetrics",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
    )


class DiffimDetectorVisitSpatiallySampledPlotsConfig(
    AnalysisBaseConfig, pipelineConnections=DiffimDetectorVisitSpatiallySampledPlotsConnections
):
    pass


class DiffimDetectorVisitSpatiallySampledPlotsTask(AnalysisPipelineTask):
    ConfigClass = DiffimDetectorVisitSpatiallySampledPlotsConfig
    _DefaultName = "DiffimDetectorVisitSpatiallySampledPlots"

    @inheritDoc(AnalysisPipelineTask)
    def parsePlotInfo(
        self, inputs: Mapping[str, Any] | None, dataId: DataCoordinate | None, connectionName: str = "data"
    ) -> Mapping[str, str]:
        """
        Notes
        -----
        This adds a 'bands' entry to `inputs`.
        """
        plotInfo = super().parsePlotInfo(inputs, dataId, connectionName=connectionName)
        plotInfo["tableName"] += f", detector: {plotInfo['detector']}"
        inputs["bands"] = plotInfo["band"]
        return plotInfo


class DiffimDetectorVisitSpatiallySampledDipoleQuiverPlotTaskConfig(
    AnalysisBaseConfig, pipelineConnections=DiffimDetectorVisitSpatiallySampledPlotsConnections
):
    pass

class DiffimDetectorVisitSpatiallySampledDipoleQuiverPlotTask(AnalysisPipelineTask):
    ConfigClass = DiffimDetectorVisitSpatiallySampledDipoleQuiverPlotTaskConfig
    _DefaultName = "DiffimDetectorVisitSpatiallySampledDipoleQuiverPlots"

    @inheritDoc(AnalysisPipelineTask)
    def parsePlotInfo(
        self, inputs: Mapping[str, Any] | None, dataId: DataCoordinate | None, connectionName: str = "data"
    ) -> Mapping[str, str]:
        """
        Notes
        -----
        This adds a 'bands' entry to `inputs`.
        """
        plotInfo = super().parsePlotInfo(inputs, dataId, connectionName=connectionName)
        plotInfo["tableName"] += f", detector: {plotInfo['detector']}"
        inputs["bands"] = plotInfo["band"]
        return plotInfo
