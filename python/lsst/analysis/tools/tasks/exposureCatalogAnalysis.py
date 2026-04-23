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

__all__ = ("ExposureCatalogAnalysisConfig", "ExposureCatalogAnalysisTask")

import numpy as np
from astropy.table import Column

from lsst.pex.config import Field, ListField
from lsst.pipe.base import connectionTypes as cT
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class ExposureCatalogAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=(),
    defaultTemplates={"inputName": ""},
):
    data = cT.Input(
        doc="Exposure Catalog to analyze.",
        name="{inputName}",
        storageClass="ExposureCatalog",
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
            storageClass=self.data.storageClass,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputTableDimensions)),
            multiple=self.data.multiple,
        )

        self.dimensions.update(frozenset(sorted(config.taskDimensions)))

        if not config.loadSkymap:
            del self.skymap
        if not config.loadCamera:
            del self.camera


class ExposureCatalogAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=ExposureCatalogAnalysisConnections
):

    inputTableDimensions = ListField[str](doc="Dimensions of the input table.", optional=False)

    taskDimensions = ListField[str](doc="Task and output dimensions.", optional=False)

    loadSkymap = Field[bool](doc="Whether to load the skymap.", default=False)

    loadCamera = Field[bool](doc="Whether to load the camera.", default=False)


class ExposureCatalogAnalysisTask(AnalysisPipelineTask):
    """An ``AnalysisPipelineTask`` that loads an ExposureCatalog passes it
    to the atools for analysis.

    It will also pass the ``camera`` and ``skyMap`` inputs to the parent run
    methods if these are requested by the config parameters.
    """

    ConfigClass = ExposureCatalogAnalysisConfig
    _DefaultName = "exposureCatalogAnalysisTask"

    def loadData(self, handle):
        # The parent class loadData does not work for ExposureCatalog.

        data = handle.get().asAstropy()
        for colname in list(data.columns.keys()):
            # An ExposureCatalog may contain elements that are themselves
            # arrays. In such cases, the arrays are expanded into multiple
            # columns named: <colname>_0, <colname>_1, etc.
            if isinstance(data[colname][0], np.ndarray):
                for index in np.arange(len(data[colname][0])):
                    values = [row[index] for row in data[colname]]
                    data.add_column(Column(name=f"{colname}_{index}", data=values))

        return data
