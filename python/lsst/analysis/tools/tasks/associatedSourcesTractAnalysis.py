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

__all__ = ("AssociatedSourcesTractAnalysisConfig", "AssociatedSourcesTractAnalysisTask")

import numpy as np
import pandas as pd
from lsst.geom import Box2D
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class AssociatedSourcesTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract", "instrument"),
    defaultTemplates={
        "outputName": "isolated_star_sources",
        "associatedSourcesInputName": "isolated_star_sources",
    },
):
    sourceCatalogs = ct.Input(
        doc="Visit based source table to load from the butler",
        name="sourceTable_visit",
        storageClass="DataFrame",
        deferLoad=True,
        dimensions=("visit", "band"),
        multiple=True,
    )

    associatedSources = ct.Input(
        doc="Table of associated sources",
        name="{associatedSourcesInputName}",
        storageClass="DataFrame",
        deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
    )

    skyMap = ct.Input(
        doc="Input definition of geometry/bbox and projection/wcs for warped exposures",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )
    camera = ct.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )


class AssociatedSourcesTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssociatedSourcesTractAnalysisConnections
):
    def setDefaults(self):
        super().setDefaults()


class AssociatedSourcesTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AssociatedSourcesTractAnalysisConfig
    _DefaultName = "associatedSourcesTractAnalysis"

    @staticmethod
    def getBoxWcs(skymap, tract):
        """Get box that defines tract boundaries."""
        tractInfo = skymap.generateTract(tract)
        wcs = tractInfo.getWcs()
        tractBox = tractInfo.getBBox()
        return tractBox, wcs

    @classmethod
    def callback(cls, inputs, dataId):
        """Callback function to be used with reconstructor."""
        return cls.prepareAssociatedSources(
            inputs["skyMap"],
            dataId["tract"],
            inputs["sourceCatalogs"],
            inputs["associatedSources"],
        )

    @classmethod
    def prepareAssociatedSources(cls, skymap, tract, sourceCatalogs, associatedSources):
        """Concatenate source catalogs and join on associated object index."""

        # Keep only sources with associations
        dataJoined = pd.concat(sourceCatalogs).merge(associatedSources, on="sourceId", how="inner")
        dataJoined.set_index("sourceId", inplace=True)

        # Determine which sources are contained in tract
        ra = np.radians(dataJoined["coord_ra"].values)
        dec = np.radians(dataJoined["coord_dec"].values)
        box, wcs = cls.getBoxWcs(skymap, tract)
        box = Box2D(box)
        x, y = wcs.skyToPixelArray(ra, dec)
        boxSelection = box.contains(x, y)

        # Keep only the sources in groups that are fully contained within the
        # tract
        dataJoined["boxSelection"] = boxSelection
        dataFiltered = dataJoined.groupby("obj_index").filter(lambda x: all(x["boxSelection"]))
        dataFiltered.drop(columns="boxSelection", inplace=True)

        return dataFiltered

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        # Load specified columns from source catalogs
        names = self.collectInputNames()
        names |= {"sourceId", "coord_ra", "coord_dec"}
        names.remove("obj_index")
        sourceCatalogs = []
        for handle in inputs["sourceCatalogs"]:
            sourceCatalogs.append(self.loadData(handle, names))
        inputs["sourceCatalogs"] = sourceCatalogs

        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId, connectionName="associatedSources")

        # TODO: make key used for object index configurable
        inputs["associatedSources"] = self.loadData(inputs["associatedSources"], ["obj_index", "sourceId"])

        data = self.callback(inputs, dataId)

        kwargs = {"data": data, "plotInfo": plotInfo, "skymap": inputs["skyMap"], "camera": inputs["camera"]}
        outputs = self.run(**kwargs)
        butlerQC.put(outputs, outputRefs)
