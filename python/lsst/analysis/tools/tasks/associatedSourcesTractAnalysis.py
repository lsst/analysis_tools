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

import lsst.pex.config as pexConfig
import numpy as np
import pandas as pd
from astropy.table import join, vstack
from lsst.drp.tasks.gbdesAstrometricFit import calculate_apparent_motion
from lsst.geom import Box2D
from lsst.pipe.base import NoWorkFound
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap
from smatch import Matcher

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
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("visit", "band"),
        multiple=True,
    )

    associatedSources = ct.Input(
        doc="Table of associated sources",
        name="{associatedSourcesInputName}",
        storageClass="ArrowAstropy",
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
    properMotionCatalogs = ct.Input(
        doc="Catalog containing proper motions.",
        name="gbdesAstrometricFit_starCatalog",
        storageClass="ArrowNumpyDict",
        dimensions=("instrument", "skymap", "tract", "physical_filter"),
        multiple=True,
        deferLoad=True,
    )
    visitTable = ct.Input(
        doc="Catalog containing visit information.",
        name="visitTable",
        storageClass="DataFrame",
        dimensions=("instrument",),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.applyProperMotions:
            self.inputs.remove("properMotionCatalogs")
            self.inputs.remove("visitTable")


class AssociatedSourcesTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssociatedSourcesTractAnalysisConnections
):
    applyProperMotions = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Apply proper motions to source positions.",
    )

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
            inputs["properMotionCatalogs"],
            inputs["visitTable"],
        )

    @classmethod
    def prepareAssociatedSources(
        cls, skymap, tract, sourceCatalogs, associatedSources, properMotionCatalogs=None, visitTable=None
    ):
        """Concatenate source catalogs and join on associated object index."""

        # Keep only sources with associations
        sourceCatalogStack = vstack(sourceCatalogs)
        dataJoined = join(sourceCatalogStack, associatedSources, keys="sourceId", join_type="inner")

        if properMotionCatalogs is not None:
            dataJoined = cls.applyProperMotionCorrections(dataJoined, properMotionCatalogs, visitTable)

        # Determine which sources are contained in tract
        ra = np.radians(dataJoined["coord_ra"])
        dec = np.radians(dataJoined["coord_dec"])
        box, wcs = cls.getBoxWcs(skymap, tract)
        box = Box2D(box)
        x, y = wcs.skyToPixelArray(ra, dec)
        boxSelection = box.contains(x, y)

        # Keep only the sources in groups that are fully contained within the
        # tract
        dataFiltered = dataJoined[boxSelection]

        return dataFiltered

    @classmethod
    def applyProperMotionCorrections(cls, dataJoined, properMotionCatalogs, visitTable):
        """Use proper motion/parallax catalogs to shift positions to median
        epoch of the visits.

        Parameters
        ----------
        dataJoined : `astropy.table.Table`
            Table containing source positions.
        properMotionCatalogs: `dict` [`pd.DataFrame`]
            Dictionary keyed by band with proper motion and parallax catalogs.
        visitTable : `pd.DataFrame`
            Table containing the MJDS of the visits.

        Returns
        -------
        dataJoined : `astropy.table.Table`
            Table containing the source positions shifted to the median visit
            epoch.
        """

        for band in np.unique(dataJoined["band"]):
            bandInd = dataJoined["band"] == band
            bandSources = dataJoined[bandInd].to_pandas()
            meanRAs = bandSources.groupby("obj_index")["coord_ra"].aggregate("mean")
            meanDecs = bandSources.groupby("obj_index")["coord_dec"].aggregate("mean")

            bandPMs = properMotionCatalogs[band]
            with Matcher(meanRAs, meanDecs) as m:
                idx, i1, i2, d = m.query_radius(
                    bandPMs["starX"], bandPMs["starY"], 0.2 / 3600, return_indices=True
                )
            import ipdb

            ipdb.set_trace()
            catRAs = np.zeros(len(meanRAs))
            catDecs = np.zeros(len(meanRAs))
            pmRAs = np.zeros(len(meanRAs))
            pmDecs = np.zeros(len(meanRAs))
            parallaxes = np.zeros(len(meanRAs))
            catRAs[i1] = bandPMs["starX"][i2]
            catDecs[i1] = bandPMs["starY"][i2]
            pmRAs[i1] = bandPMs["starPMx"][i2]
            pmDecs[i1] = bandPMs["starPMy"][i2]
            parallaxes[i1] = bandPMs["starParallax"][i2]
            pmDf = pd.DataFrame(
                {
                    "ra": catRAs,
                    "dec": catDecs,
                    "pmRA": pmRAs,
                    "pmDec": pmDecs,
                    "parallax": parallaxes,
                    "obj_index": meanRAs.index,
                }
            )

            dataWithPM = pd.merge(bandSources, pmDf, on="obj_index", how="left")

            visits = bandSources["visit"].unique()
            mjds = [visitTable.loc[visit]["expMidptMJD"] for visit in visits]
            mjdDf = pd.DataFrame({"MJD": mjds, "visit": visits})
            dataWithMJD = pd.merge(dataWithPM, mjdDf, on="visit", how="left")
            medianMJD = np.median(mjds)
            raCorrection, decCorrection = calculate_apparent_motion(dataWithMJD, medianMJD)

            dataJoined["coord_ra"][bandInd] = dataWithMJD["coord_ra"] - raCorrection.value
            dataJoined["coord_dec"][bandInd] = dataWithMJD["coord_dec"] - decCorrection.value

        return dataJoined

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

        if self.config.applyProperMotions:
            properMotions = {}
            for pmCatRef in inputs["properMotionCatalogs"]:
                pmCat = pmCatRef.get(
                    parameters={"columns": ["starX", "starY", "starPMx", "starPMy", "starParallax"]}
                )
                properMotions[pmCatRef.dataId["band"]] = pd.DataFrame(pmCat)
            inputs["properMotionCatalogs"] = properMotions
        else:
            inputs["properMotionCatalogs"] = None
            inputs["visitTable"] = None

        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId, connectionName="associatedSources")

        # TODO: make key used for object index configurable
        inputs["associatedSources"] = self.loadData(inputs["associatedSources"], ["obj_index", "sourceId"])

        if len(inputs["associatedSources"]) == 0:
            raise NoWorkFound(f"No associated sources in tract {dataId.tract.id}")

        data = self.callback(inputs, dataId)

        kwargs = {"data": data, "plotInfo": plotInfo, "skymap": inputs["skyMap"], "camera": inputs["camera"]}
        outputs = self.run(**kwargs)
        butlerQC.put(outputs, outputRefs)
