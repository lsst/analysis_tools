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

__all__ = (
    "WholeTractImageAnalysisConfig",
    "WholeTractImageAnalysisTask",
)

from typing import Any, Mapping

from lsst.daf.butler import DataCoordinate
from lsst.pipe.base import (
    InputQuantizedConnection,
    OutputQuantizedConnection,
    QuantumContext,
)
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class WholeTractImageAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract", "band"),
    defaultTemplates={
        "inputName": "deep_coaddBin",
        "outputName": "deepTract_PostageStamp",
    },
):
    data = ct.Input(
        doc="Binned deepCoadd calibrated exposures to read from the butler.",
        name="{inputName}",
        storageClass="ExposureF",
        deferLoad=True,
        dimensions=(
            "skymap",
            "tract",
            "patch",
            "band",
        ),
        multiple=True,
    )

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )


class WholeTractImageAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=WholeTractImageAnalysisConnections
):
    pass


class WholeTractImageAnalysisTask(AnalysisPipelineTask):

    ConfigClass = WholeTractImageAnalysisConfig
    _DefaultName = "wholeTractImageAnalysisTask"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)

        try:
            inputData = inputs.pop("data")
        except KeyError:
            raise RuntimeError("'data' is a required input connection, but is not defined.")

        inputNames = {"mask"}
        inputNames.update(self.collectInputNames())

        keyedData = dict()
        for inputName in inputNames:
            keyedData[inputName] = dict()
            for handle in inputData:
                keyedData[inputName][handle.dataId["patch"]] = handle.get(component=inputName)

        outputs = self.run(
            data=keyedData,
            plotInfo=plotInfo,
            tractId=dataId["tract"],
            skymap=inputs["skymap"],
            bands=dataId["band"],
        )

        self.putByBand(butlerQC, outputs, outputRefs)

    def parsePlotInfo(
        self, inputs: Mapping[str, Any] | None, dataId: DataCoordinate | None, connectionName: str = "data"
    ) -> Mapping[str, str]:
        """Parse the inputs and dataId to get the information needed to
        to add to the figure. The parent class parsePlotInfo cannot be
        used becuase it assumes a single input dataset, as opposed to the
        multiple datasets used by this analysis task.

        Parameters
        ----------
        inputs: `dict`
            The inputs to the task
        dataCoordinate: `lsst.daf.butler.DataCoordinate`
            The dataId that the task is being run on.
        connectionName: `str`, optional
            Name of the input connection to use for determining table name.

        Returns
        -------
        plotInfo : `dict`
        """

        if inputs is None:
            tableName = ""
            run = ""
        else:
            tableName = inputs[connectionName][0].ref.datasetType.name
            run = inputs[connectionName][0].ref.run

        # Initialize the plot info dictionary
        plotInfo = {"tableName": tableName, "run": run}

        self._populatePlotInfoWithDataId(plotInfo, dataId)
        return plotInfo
