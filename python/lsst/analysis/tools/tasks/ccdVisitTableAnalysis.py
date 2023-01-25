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

from lsst.cp.pipe._lookupStaticCalibration import lookupStaticCalibration
from lsst.pex.config import Field
from lsst.pipe.base import connectionTypes as cT
from lsst.skymap import BaseSkyMap

from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CcdVisitTableAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "instrument"),
    defaultTemplates={"outputName": "ccdVisitTable"},
):
    data = cT.Input(
        doc="Collection based source table to load from the butler.",
        name="ccdVisitTable",
        storageClass="DataFrame",
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
        lookupFunction=lookupStaticCalibration,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        # No metrics are computed for this task, so remove output dataset.
        self.outputs.remove("metrics")
        if not config.introspectMakeWarpConfig:
            self.inputs.remove("makeWarpConfig")


class CcdVisitTableAnalysisConfig(AnalysisBaseConfig, pipelineConnections=CcdVisitTableAnalysisConnections):
    introspectMakeWarpConfig = Field[bool](
        doc="Whether to introspect the makeWarp_config dataset to obtain the actual "
        "maxEllipResidual and maxScalesSizeScatter thresholds in this run?  Set to "
        "False if makeWarp has not yet been run on given collection.",
        default=True,
    )


class CcdVisitTableAnalysisTask(AnalysisPipelineTask):
    ConfigClass = CcdVisitTableAnalysisConfig
    _DefaultName = "ccdVisitTableAnalysisTask"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class.
        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)
        data = self.loadData(inputs["data"])
        skymap = None if "skymap" not in inputs.keys() else inputs["skymap"]
        camera = None if "camera" not in inputs.keys() else inputs["camera"]
        calibrateConfig = inputs["calibrateConfig"]
        makeWarpConfig = None if "makeWarpConfig" not in inputs.keys() else inputs["makeWarpConfig"]

        outputs = self.run(
            data=data,
            plotInfo=plotInfo,
            camera=camera,
            skymap=skymap,
            calibrateConfig=calibrateConfig,
            makeWarpConfig=makeWarpConfig,
        )
        butlerQC.put(outputs, outputRefs)
