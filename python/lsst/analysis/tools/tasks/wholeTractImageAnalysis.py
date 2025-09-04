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
    "MakeBinnedCoaddConfig",
    "MakeBinnedCoaddTask",
)

from typing import Any, Mapping

import lsst.pipe.base as pipeBase
from lsst.daf.butler import DataCoordinate
from lsst.ip.isr.binImageDataTask import binImageData
from lsst.pex.config import Field
from lsst.pipe.base import (
    InputQuantizedConnection,
    OutputQuantizedConnection,
    PipelineTask,
    PipelineTaskConfig,
    PipelineTaskConnections,
    QuantumContext,
)
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class WholeTractImageAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract", "band"),
    defaultTemplates={
        "coaddName": "deep",
    },
):
    data = ct.Input(
        doc="Binned coadd image data to read from the butler.",
        name="{coaddName}Coadd_calexp_bin",
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

    def __init__(self, *, config=None):
        """Customize the storageClass for a specific instance. This enables it
        to be dynamically set at runtime, allowing the task to work with
        different types of image-like data.

        Parameters
        ----------
        config : `WholeTractImageAnalysisConfig`
            A config for `WholeTractImageAnalysisTask`.
        """
        super().__init__(config=config)
        if config and config.dataStorageClass != self.data.storageClass:
            self.data = ct.Input(
                name=self.data.name,
                doc=self.data.doc,
                storageClass=config.dataStorageClass,
                dimensions=self.data.dimensions,
                deferLoad=self.data.deferLoad,
                multiple=self.data.multiple,
            )


class WholeTractImageAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=WholeTractImageAnalysisConnections
):
    dataStorageClass = Field(
        default="ExposureF",
        dtype=str,
        doc=(
            "Override the storageClass of the input data. "
            "Must be of type `Image`, `MaskedImage` or `Exposure`, or one of their subtypes."
        ),
    )


class WholeTractImageAnalysisTask(AnalysisPipelineTask):

    ConfigClass = WholeTractImageAnalysisConfig
    _DefaultName = "wholeTractImageAnalysis"

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

        keyedData = dict()
        if "Exposure" in self.config.dataStorageClass:
            inputNames = {"mask"}
            inputNames.update(self.collectInputNames())
            for inputName in inputNames:
                keyedData[inputName] = dict()
                for handle in inputData:
                    keyedData[inputName][handle.dataId["patch"]] = handle.get(component=inputName)
        elif "Image" in self.config.dataStorageClass:
            keyedData["image"] = dict()
            for handle in inputData:
                image = handle.get()
                keyedData["image"][handle.dataId["patch"]] = image
        else:
            raise TypeError("'data' must be of type Image, MaskedImage, Exposure, or one of their subtypes")

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


class MakeBinnedCoaddConnections(
    PipelineTaskConnections,
    dimensions=("skymap", "tract", "patch", "band"),
    defaultTemplates={"coaddName": "deep"},
):

    coadd = ct.Input(
        doc="Input coadd image data to bin.",
        name="{coaddName}Coadd_calexp",
        storageClass="ExposureF",
        dimensions=("skymap", "tract", "patch", "band"),
        deferLoad=True,
    )
    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )
    binnedCoadd = ct.Output(
        doc="Binned coadd image data.",
        name="{coaddName}Coadd_calexp_bin",
        storageClass="ExposureF",
        dimensions=("skymap", "tract", "patch", "band"),
    )

    def __init__(self, *, config=None):
        """Customize the storageClass for a specific instance.
        This enables it to be dynamically set at runtime, allowing
        the task to work with different types of image-like data.

        Parameters
        ----------
        config : `MakeBinnedCoaddConfig`
            A config for `MakeBinnedCoaddTask`.
        """
        super().__init__(config=config)
        if config and config.coaddStorageClass != self.coadd.storageClass:
            self.coadd = ct.Input(
                name=self.coadd.name,
                doc=self.coadd.doc,
                storageClass=config.coaddStorageClass,
                dimensions=self.coadd.dimensions,
                deferLoad=self.coadd.deferLoad,
            )
            self.binnedCoadd = ct.Output(
                name=self.binnedCoadd.name,
                doc=self.binnedCoadd.doc,
                storageClass=config.coaddStorageClass,
                dimensions=self.binnedCoadd.dimensions,
            )


class MakeBinnedCoaddConfig(PipelineTaskConfig, pipelineConnections=MakeBinnedCoaddConnections):
    """Config for MakeBinnedCoaddTask"""

    doBinInnerBBox = Field[bool](
        doc=(
            "Retrieve and bin the coadd image data within the patch Inner Bounding Box, ",
            "thereby excluding the regions that overlap neighboring patches.",
        ),
        default=False,
    )
    binFactor = Field[int](
        doc="Binning factor applied to both spatial dimensions.",
        default=8,
        check=lambda x: x > 1,
    )
    coaddStorageClass = Field(
        default="ExposureF",
        dtype=str,
        doc=(
            "Override the storageClass of the input and binned coadd image data. "
            "Must be of type `Image`, `MaskedImage`, or `Exposure`, or one of their subtypes."
        ),
    )


class MakeBinnedCoaddTask(PipelineTask):

    ConfigClass = MakeBinnedCoaddConfig
    _DefaultName = "makeBinnedCoadd"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        """Takes coadd image data and bins it by the factor specified in
        self.config.binFactor. This task uses the binImageData function
        defined in ip_isr, but adds the option to only retrieve and bin the
        data contained within the patch's inner bounding box.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
            A butler which is specialized to operate in the context of a
            `lsst.daf.butler.Quantum`.
        inputRefs : `lsst.pipe.base.InputQuantizedConnection`
            Data structure containing named attributes 'coadd' and 'skymap'.
            The values of these attributes are the corresponding
            `lsst.daf.butler.DatasetRef` objects defined in the corresponding
            `PipelineTaskConnections` class.
        outputRefs : `lsst.pipe.base.OutputQuantizedConnection`
            Datastructure containing named attribute 'binnedCoadd'.
            The value of this attribute is the corresponding
            `lsst.daf.butler.DatasetRef` object defined in the corresponding
            `PipelineTaskConnections` class.
        """

        inputs = butlerQC.get(inputRefs)
        coaddRef = inputs["coadd"]

        if self.config.doBinInnerBBox:
            skymap = inputs["skymap"]
            tractId = butlerQC.quantum.dataId["tract"]
            patchId = butlerQC.quantum.dataId["patch"]
            tractInfo = skymap.generateTract(tractId)
            bbox = tractInfo.getPatchInfo(patchId).getInnerBBox()

            coadd = coaddRef.get(parameters={"bbox": bbox})
        else:
            coadd = coaddRef.get()

        binnedCoadd = binImageData(coadd, self.config.binFactor)

        butlerQC.put(pipeBase.Struct(binnedCoadd=binnedCoadd), outputRefs)
