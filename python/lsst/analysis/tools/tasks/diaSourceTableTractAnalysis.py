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

__all__ = ("DiaSourceTableTractAnalysisConfig", "DiaSourceTableTractAnalysisTask")

from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask
from ..atools import diaSourceTableTractMetrics


class DiaSourceTableTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={
        "consolidateAssocDiaSourceTableInputName": "diaSourceTable_tract",
    },
):
    data = ct.Input(
        doc="Table of per-tract DiaSources to load from the butler",
        name="{consolidateAssocDiaSourceTableInputName}",
        storageClass="DataFrame",
        deferLoad=True,
        dimensions=(
            "skymap",
            "tract",
        ),
    )

    skyMap = ct.Input(
        doc="Input definition of geometry/bbox and projection/wcs",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )


class DiaSourceTableTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiaSourceTableTractAnalysisConnections
):
    def setDefaults(self):
        """All the metrics are turned on in this task by default,
        so it is not necessary to configure each one in a pipeline
        """
        super().setDefaults()
        self.atools.NumDiaSources = diaSourceTableTractMetrics.NumDiaSourcesMetric()
        self.atools.NumStreakDiaSources = diaSourceTableTractMetrics.NumStreakDiaSourcesMetric()
        self.atools.NumStreakCenterDiaSources = diaSourceTableTractMetrics.NumStreakCenterDiaSourcesMetric()


class DiaSourceTableTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiaSourceTableTractAnalysisConfig
    _DefaultName = "DiaSourceTableTractAnalysis"
