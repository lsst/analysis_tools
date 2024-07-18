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

import astropy.time
import astropy.units as u
import lsst.pex.config as pexConfig
import numpy as np
import pandas as pd
from astropy.table import Table, join, vstack
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
    astrometricCorrectionCatalogs = ct.Input(
        doc="Catalog containing proper motion and parallax.",
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

        if not config.applyAstrometricCorrections:
            self.inputs.remove("astrometricCorrectionCatalogs")
            self.inputs.remove("visitTable")


class AssociatedSourcesTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssociatedSourcesTractAnalysisConnections
):
    applyAstrometricCorrections = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Apply proper motion and parallax corrections to source positions.",
    )
    matchingRadius = pexConfig.Field(
        dtype=float,
        default=0.2,
        doc=(
            "Radius in mas with which to match the mean positions of the sources with the positions in the"
            " astrometricCorrectionCatalogs."
        ),
    )
    astrometricCorrectionParameters = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        # TODO: DM-45845 Update default names when catalog gets updated.
        default={
            "ra": "starX",
            "dec": "starY",
            "pmRA": "starPMx",
            "pmDec": "starPMy",
            "parallax": "starParallax",
        },
        doc="Column names for position and motion parameters in the astrometric correction catalogs.",
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

    def callback(self, inputs, dataId):
        """Callback function to be used with reconstructor."""
        return self.prepareAssociatedSources(
            inputs["skyMap"],
            dataId["tract"],
            inputs["sourceCatalogs"],
            inputs["associatedSources"],
            inputs["astrometricCorrectionCatalogs"],
            inputs["visitTable"],
        )

    def prepareAssociatedSources(
        self,
        skymap,
        tract,
        sourceCatalogs,
        associatedSources,
        astrometricCorrectionCatalogs=None,
        visitTable=None,
    ):
        """Concatenate source catalogs and join on associated object index."""

        # Keep only sources with associations
        sourceCatalogStack = vstack(sourceCatalogs)
        dataJoined = join(sourceCatalogStack, associatedSources, keys="sourceId", join_type="inner")

        if astrometricCorrectionCatalogs is not None:
            self.applyAstrometricCorrections(dataJoined, astrometricCorrectionCatalogs, visitTable)

        # Determine which sources are contained in tract
        ra = np.radians(dataJoined["coord_ra"])
        dec = np.radians(dataJoined["coord_dec"])
        box, wcs = self.getBoxWcs(skymap, tract)
        box = Box2D(box)
        x, y = wcs.skyToPixelArray(ra, dec)
        boxSelection = box.contains(x, y)

        # Keep only the sources in groups that are fully contained within the
        # tract
        dataFiltered = dataJoined[boxSelection]

        return dataFiltered

    def applyAstrometricCorrections(self, dataJoined, astrometricCorrectionCatalogs, visitTable):
        """Use proper motion/parallax catalogs to shift positions to median
        epoch of the visits.

        Parameters
        ----------
        dataJoined : `astropy.table.Table`
            Table containing source positions, which will be modified in place.
        astrometricCorrectionCatalogs: `dict` [`pd.DataFrame`]
            Dictionary keyed by band with proper motion and parallax catalogs.
        visitTable : `pd.DataFrame`
            Table containing the MJDs of the visits.
        """
        for band in np.unique(dataJoined["band"]):
            bandInd = dataJoined["band"] == band
            bandSources = dataJoined[bandInd]
            # Add key for sorting below.
            bandSources["__index__"] = np.arange(len(bandSources))
            bandSourcesDf = bandSources.to_pandas()
            meanRAs = bandSourcesDf.groupby("obj_index")["coord_ra"].aggregate("mean")
            meanDecs = bandSourcesDf.groupby("obj_index")["coord_dec"].aggregate("mean")

            bandPMs = astrometricCorrectionCatalogs[band]
            with Matcher(meanRAs, meanDecs) as m:
                idx, i1, i2, d = m.query_radius(
                    bandPMs[self.config.astrometricCorrectionParameters["ra"]],
                    bandPMs[self.config.astrometricCorrectionParameters["dec"]],
                    (self.config.matchingRadius * u.mas).to(u.degree),
                    return_indices=True,
                )

            catRAs = np.zeros_like(meanRAs)
            catDecs = np.zeros_like(meanRAs)
            pmRAs = np.zeros_like(meanRAs)
            pmDecs = np.zeros_like(meanRAs)
            parallaxes = np.zeros(len(meanRAs))
            catRAs[i1] = bandPMs[self.config.astrometricCorrectionParameters["ra"]][i2]
            catDecs[i1] = bandPMs[self.config.astrometricCorrectionParameters["dec"]][i2]
            pmRAs[i1] = bandPMs[self.config.astrometricCorrectionParameters["pmRA"]][i2]
            pmDecs[i1] = bandPMs[self.config.astrometricCorrectionParameters["pmDec"]][i2]
            parallaxes[i1] = bandPMs[self.config.astrometricCorrectionParameters["parallax"]][i2]

            pmDf = Table(
                {
                    "ra": catRAs * u.degree,
                    "dec": catDecs * u.degree,
                    "pmRA": pmRAs * u.mas / u.yr,
                    "pmDec": pmDecs * u.mas / u.yr,
                    "parallax": parallaxes * u.mas,
                    "obj_index": meanRAs.index,
                }
            )

            dataWithPM = join(bandSources, pmDf, keys="obj_index", join_type="left")

            visits = bandSourcesDf["visit"].unique()
            mjds = [visitTable.loc[visit]["expMidptMJD"] for visit in visits]
            mjdTable = Table(
                [astropy.time.Time(mjds, format="mjd", scale="tai"), visits], names=["MJD", "visit"]
            )
            dataWithMJD = join(dataWithPM, mjdTable, keys="visit", join_type="left")
            # After astropy 7.0, it should be possible to use "keep_order=True"
            # in the join and avoid sorting.
            dataWithMJD.sort("__index__")
            medianMJD = astropy.time.Time(np.median(mjds), format="mjd", scale="tai")

            raCorrection, decCorrection = calculate_apparent_motion(dataWithMJD, medianMJD)

            dataJoined["coord_ra"][bandInd] = dataWithMJD["coord_ra"] - raCorrection.value
            dataJoined["coord_dec"][bandInd] = dataWithMJD["coord_dec"] - decCorrection.value

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

        if self.config.applyAstrometricCorrections:
            astrometricCorrections = {}
            for pmCatRef in inputs["astrometricCorrectionCatalogs"]:
                pmCat = pmCatRef.get(
                    parameters={"columns": self.config.astrometricCorrectionParameters.values()}
                )
                astrometricCorrections[pmCatRef.dataId["band"]] = pd.DataFrame(pmCat)
            inputs["astrometricCorrectionCatalogs"] = astrometricCorrections
        else:
            inputs["astrometricCorrectionCatalogs"] = None
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
