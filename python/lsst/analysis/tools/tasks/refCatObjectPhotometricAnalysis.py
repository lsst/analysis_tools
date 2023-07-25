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

__all__ = ("RefCatObjectPhotometricAnalysisConfig", "RefCatObjectPhotometricAnalysisTask")

from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class RefCatObjectPhotometricAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={"outputName": "objectTable_tract_ps1_pv3_3pi_20170110_match"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="objectTable_tract_ps1_pv3_3pi_20170110_match",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("skymap", "tract"),
    )

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )


class RefCatObjectPhotometricAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=RefCatObjectPhotometricAnalysisConnections
):
    pass


class RefCatObjectPhotometricAnalysisTask(AnalysisPipelineTask):
    """Make plots and metrics using a table of objects
    matched to photometric reference catalog sources.
    These plots compare the photometry in each of the
    specified bands.
    """

    ConfigClass = RefCatObjectPhotometricAnalysisConfig
    _DefaultName = "refCatObjectPhotometricAnalysisTask"
