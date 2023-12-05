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

__all__ = ("CcdVisitTableAnalysisConfig", "CcdVisitTableAnalysisTask")

from lsst.pex.config import Field
from lsst.pipe.base import connectionTypes as cT
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CcdVisitTableAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "instrument"),
    defaultTemplates={"inputName": "ccdVisitTable", "outputName": "ccdVisitTable"},
):
    data = cT.Input(
        doc="Collection based source table to load from the butler.",
        name="{inputName}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument",),
    )
    calibrateConfig = cT.Input(
        doc="Collection based calibrate task config to load from the butler.",
        name="calibrate_config",
        storageClass="Config",
        dimensions=(),
    )
    makeWarpConfig = cT.Input(
        doc="Collection based makeWarp task config to load from the butler.",
        name="makeWarp_config",
        storageClass="Config",
        dimensions=(),
    )
    skymap = cT.Input(
        doc="Input definition of geometry/bbox and projection/wcs for warped exposures.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )
    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        if not config.introspectMakeWarpConfig:
            self.inputs.remove("makeWarpConfig")


class CcdVisitTableAnalysisConfig(AnalysisBaseConfig, pipelineConnections=CcdVisitTableAnalysisConnections):
    introspectMakeWarpConfig = Field[bool](
        doc="Whether to introspect the makeWarp_config dataset to obtain the actual "
        "maxEllipResidual and maxScalesSizeScatter thresholds in this run?  Set to "
        "True only if makeWarp has been run on given collection.  When False, a default "
        "set of threshold values (meant to reflect current stack defaults) will be set "
        "for reference.",
        default=False,
    )


class CcdVisitTableAnalysisTask(AnalysisPipelineTask):
    ConfigClass = CcdVisitTableAnalysisConfig
    _DefaultName = "ccdVisitTableAnalysisTask"
