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

__all__ = ("GenericTableAnalysisConfig", "GenericTableAnalysisTask")

from lsst.pex.config import Field, ListField
from lsst.pipe.base import connectionTypes as cT
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class GenericTableAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=(),
    defaultTemplates={"inputName": ""},
):
    data = cT.Input(
        doc="Table to analyze.",
        name="{inputName}",
        storageClass="ArrowAstropy",
        dimensions=(),
        deferLoad=True,
    )

    skymap = cT.Input(
        doc="The skymap that covers the tract that the data is from.",
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

        self.data = cT.Input(
            doc=self.data.doc,
            name=self.data.name,
            storageClass=config.inputTableStorageClass,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputTableDimensions)),
            multiple=self.data.multiple,
        )

        self.dimensions.update(frozenset(sorted(config.taskDimensions)))

        if not config.loadSkymap:
            del self.skymap
        if not config.loadCamera:
            del self.camera


class GenericTableAnalysisConfig(AnalysisBaseConfig, pipelineConnections=GenericTableAnalysisConnections):

    inputTableStorageClass = Field[str](doc="Storage class of input table.", optional=False)

    inputTableDimensions = ListField[str](doc="Dimensions of the input table.", optional=False)

    taskDimensions = ListField[str](doc="Task dimensions", optional=False)

    loadSkymap = Field[bool](doc="Whether to load the skymap", default=False, optional=True)

    loadCamera = Field[bool](doc="Whether to load the camera", default=False, optional=True)


class GenericTableAnalysisTask(AnalysisPipelineTask):
    """An ``AnalysisPipelineTask`` that loads a generic table and passes it
    to the atools for analysis.

    Supports ``ExposureCatalog`` and ``ArrowAstropy`` storage classes. It
    will also pass the ``camera`` and ``skyMap`` inputs to the parent run
    methods if these are requested by the config parameters.
    """

    ConfigClass = GenericTableAnalysisConfig
    _DefaultName = "genericTableAnalysisTask"

    def loadData(self, handle, names=None):

        storageClassName = handle.ref.datasetType.storageClass_name
        if storageClassName == "ArrowAstropy":
            data = super().loadData(handle)
        elif storageClassName == "ExposureCatalog":
            data = handle.get()
        else:
            raise ValueError(f"Unsupported input table storage class: {storageClassName}")

        return data
