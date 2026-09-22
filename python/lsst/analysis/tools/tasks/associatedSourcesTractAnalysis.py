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

__all__ = (
    "AssociatedSourcesTractAnalysisConfig",
    "AssociatedSourcesTractAnalysisTask",
    "AssociatedSourcesHealpix3AnalysisTask",
)

import astropy.time
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, vstack
from scipy.spatial import KDTree

import lsst.pex.config as pexConfig
from lsst.daf.butler import DatasetProvenance
from lsst.drp.tasks.gbdesAstrometricFit import calculate_apparent_motion
from lsst.pipe.base import NoWorkFound
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap
from lsst.sphgeom import HealpixPixelization

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class AssociatedSourcesTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract", "instrument"),
    defaultTemplates={
        "outputName": "isolated_star_presources",
        "associatedSourcesInputName": "isolated_star_presources",
        "associatedSourceIdsInputName": "isolated_star_presource_associations",
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

    associatedSourceIds = ct.Input(
        doc="Table containing unique ids for the associated sources",
        name="{associatedSourceIdsInputName}",
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
    astrometricCorrectionCatalog = ct.Input(
        doc="Catalog with proper motion and parallax information.",
        name="isolated_star_stellar_motions",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
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
            self.inputs.remove("astrometricCorrectionCatalog")
            self.inputs.remove("visitTable")


class AssociatedSourcesTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssociatedSourcesTractAnalysisConnections
):
    applyAstrometricCorrections = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Apply proper motion and parallax corrections to source positions.",
    )
    astrometricCorrectionParameters = pexConfig.DictField(
        keytype=str,
        itemtype=str,
        default={
            "ra": "ra",
            "dec": "dec",
            "pmRA": "raPM",
            "pmDec": "decPM",
            "parallax": "parallax",
            "isolated_star_id": "isolated_star_id",
        },
        doc="Column names for position and motion parameters in the astrometric correction catalogs.",
    )
    maxVisitCount = pexConfig.Field(
        dtype=int,
        default=None,
        doc="Maximum number of visits to use in calculating the metrics.",
        optional=True,
    )
    maxObjects = pexConfig.Field(
        dtype=int,
        default=None,
        doc="Maximum number of associated objects to use in calculating the metrics.",
        optional=True,
    )


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
            inputs["sourceCatalogs"],
            inputs["associatedSources"],
            inputs["associatedSourceIds"],
            inputs["astrometricCorrectionCatalog"],
            inputs["visitTable"],
        )

    def prepareAssociatedSources(
        self,
        sourceCatalogs,
        associatedSources,
        associatedSourceIds,
        astrometricCorrectionCatalog=None,
        visitTable=None,
    ):
        """Concatenate source catalogs and join on associated source IDs."""
        rng = np.random.default_rng()

        # Strip any provenance from tables before merging to prevent
        # warnings from conflicts being issued by astropy.utils.merge.
        DatasetProvenance.strip_provenance_from_flat_dict(associatedSources.meta)
        DatasetProvenance.strip_provenance_from_flat_dict(associatedSourceIds.meta)

        # associatedSource["obj_index"] refers to the corresponding index (row)
        # in associatedSourceIds.
        index = associatedSources["obj_index"]
        associatedSources["isolated_star_id"] = associatedSourceIds["isolated_star_id"][index]

        if self.config.maxObjects:
            objectChoice = rng.permutation(associatedSourceIds["isolated_star_id"])[: self.config.maxObjects]
            objectChoice.sort()
            sub1 = np.clip(
                np.searchsorted(objectChoice, associatedSources["isolated_star_id"]),
                0,
                len(objectChoice) - 1,
            )
            matched = objectChoice[sub1] == associatedSources["isolated_star_id"]
            associatedSources = associatedSources[matched]

        trimmedSourceCatalogs = []
        fullCatLen = 0
        # It would be preferable to use astropy's built in functions
        # but they are too slow so instead we use a numpy searchsorted
        # maneuver.
        sortedAssocSources = associatedSources["sourceId"].copy()
        assocSourcesSort = associatedSources["sourceId"].argsort()
        sortedAssocSources.sort()
        nAssocSources = len(sortedAssocSources)
        colsNeeded = list(self.collectInputNames())
        # Only get the columns needed for the source catalogues.
        # The isolated_star_id and the obj_index are added later
        # from other tables so remove these from the list. Also
        # add the coord_ra and coord_dec as well because this bit
        # of code needs it even if it isn't requested by a
        # downstream atool.
        if "isolated_star_id" in colsNeeded:
            colsNeeded.remove("isolated_star_id")
        if "obj_index" in colsNeeded:
            colsNeeded.remove("obj_index")
        colsNeeded += ["sourceId", "coord_ra", "coord_dec"]

        if self.config.maxVisitCount:
            sourceCatalogs = rng.permutation(sourceCatalogs)[: self.config.maxVisitCount]

        for sourceCatalogRef in sourceCatalogs:
            sourceCatalog = sourceCatalogRef.get(parameters={"columns": set(colsNeeded)})
            DatasetProvenance.strip_provenance_from_flat_dict(sourceCatalog.meta)

            sub = np.clip(
                np.searchsorted(sortedAssocSources, sourceCatalog["sourceId"]), 0, nAssocSources - 1
            )
            sourceCatalogInds = sortedAssocSources[sub] == sourceCatalog["sourceId"]
            assocCatalogInds = sub[sourceCatalogInds]

            # Keep only the sources in groups that are fully contained within
            # the tract by matching to the associated sources table
            trimmedSourceCatalogs.append(
                hstack(
                    [associatedSources[assocSourcesSort][assocCatalogInds], sourceCatalog[sourceCatalogInds]]
                )
            )
            fullCatLen += np.sum(sourceCatalogInds)

        columns = trimmedSourceCatalogs[0].columns
        dtypes = trimmedSourceCatalogs[0].dtype
        zeros = np.zeros((fullCatLen, len(columns)))
        fullCat = Table(data=zeros, names=columns, dtype=dtypes)
        n = 0
        for trimmedSourceCatalog in trimmedSourceCatalogs:
            fullCat[n : n + len(trimmedSourceCatalog)] = trimmedSourceCatalog
            n += len(trimmedSourceCatalog)

        if (astrometricCorrectionCatalog is not None) and (len(fullCat) != 0):
            self.applyAstrometricCorrections(fullCat, astrometricCorrectionCatalog, visitTable)

        # Keep only finite ras and decs
        keep = np.isfinite(fullCat["coord_ra"]) & np.isfinite(fullCat["coord_dec"])
        return fullCat[keep]

    def applyAstrometricCorrections(self, dataJoined, astrometricCorrectionCatalog, visitTable):
        """Use proper motion/parallax catalogs to shift positions to median
        epoch of the visits.

        Parameters
        ----------
        dataJoined : `astropy.table.Table`
            Table containing source positions, which will be modified in place.
        astrometricCorrectionCatalog : `astropy.table.Table`
            Proper motion and parallax catalog.
        visitTable : `pd.DataFrame`
            Table containing the MJDs of the visits.
        """
        if visitTable.index.name is None:
            # The expected index may or may not be set, depending on whether
            # the table was written originally as a DataFrame or something else
            # Parquet-friendly.
            visitTable.set_index("visitId", inplace=True)

        # Get the stellar motion catalog into the right format:
        for key, value in self.config.astrometricCorrectionParameters.items():
            astrometricCorrectionCatalog.rename_column(value, key)
        astrometricCorrectionCatalog["ra"] *= u.degree
        astrometricCorrectionCatalog["dec"] *= u.degree
        astrometricCorrectionCatalog["pmRA"] *= u.mas / u.yr
        astrometricCorrectionCatalog["pmDec"] *= u.mas / u.yr
        astrometricCorrectionCatalog["parallax"] *= u.mas

        # Again using astropy join would have been great but this is four
        # times faster
        lenAstroCorrCat = len(astrometricCorrectionCatalog)
        tree = KDTree(astrometricCorrectionCatalog["isolated_star_id"].reshape(lenAstroCorrCat, 1))
        _, inds = tree.query(
            dataJoined["isolated_star_id"].reshape(len(dataJoined), 1), distance_upper_bound=0.5
        )
        ids = inds < lenAstroCorrCat

        dataWithPM = hstack([dataJoined[ids], astrometricCorrectionCatalog[inds[ids]]])

        mjds = visitTable.loc[dataWithPM["visit"]]["expMidptMJD"]
        times = astropy.time.Time(mjds, format="mjd", scale="tai")
        dataWithPM["MJD"] = times
        medianMJD = astropy.time.Time(np.median(mjds), format="mjd", scale="tai")

        raCorrection, decCorrection = calculate_apparent_motion(dataWithPM, medianMJD)

        dataJoined["coord_ra"][ids] = dataWithPM["coord_ra"] - raCorrection.value
        dataJoined["coord_dec"][ids] = dataWithPM["coord_dec"] - decCorrection.value

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        # Load specified columns from source catalogs
        names = self.collectInputNames()
        names |= {"sourceId", "coord_ra", "coord_dec"}
        for item in ["obj_index", "isolated_star_id"]:
            if item in names:
                names.remove(item)

        if self.config.applyAstrometricCorrections:
            astrometricCorrections = inputs["astrometricCorrectionCatalog"].get(
                parameters={"columns": self.config.astrometricCorrectionParameters.values()}
            )
            inputs["astrometricCorrectionCatalog"] = astrometricCorrections
        else:
            inputs["astrometricCorrectionCatalog"] = None
            inputs["visitTable"] = None

        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId, connectionName="associatedSources")

        # TODO: make key used for object index configurable
        inputs["associatedSources"] = self.loadData(inputs["associatedSources"], ["obj_index", "sourceId"])
        inputs["associatedSourceIds"] = self.loadData(inputs["associatedSourceIds"], ["isolated_star_id"])

        if len(inputs["associatedSources"]) == 0:
            raise NoWorkFound(f"No associated sources in tract {dataId.tract.id}")

        data = self.callback(inputs, dataId)

        kwargs = {"data": data, "plotInfo": plotInfo, "skymap": inputs["skyMap"], "camera": inputs["camera"]}
        outputs = self.run(**kwargs)
        self.putByBand(butlerQC, outputs, outputRefs)


class AssociatedSourcesHealpix3AnalysisConnections(
    AssociatedSourcesTractAnalysisConnections,
    dimensions=("healpix3", "instrument"),
):
    associatedSources = ct.Input(
        doc="Table of associated sources",
        name="{associatedSourcesInputName}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
        multiple=True,
    )

    associatedSourceIds = ct.Input(
        doc="Table containing unique ids for the associated sources",
        name="{associatedSourceIdsInputName}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
        multiple=True,
    )
    astrometricCorrectionCatalog = ct.Input(
        doc="Catalog with proper motion and parallax information.",
        name="isolated_star_stellar_motions",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
        multiple=True,
    )


class AssociatedSourcesHealpix3AnalysisConfig(
    AssociatedSourcesTractAnalysisConfig, pipelineConnections=AssociatedSourcesHealpix3AnalysisConnections
):
    pass


class AssociatedSourcesHealpix3AnalysisTask(AssociatedSourcesTractAnalysisTask):
    ConfigClass = AssociatedSourcesHealpix3AnalysisConfig
    _DefaultName = "associatedSourcesHealpix3Analysis"

    def getHealpixOverlap(self, sources, sourceIds, pixelId, astrometricCorrections=None):

        pixelization = HealpixPixelization(3)
        pixelRegion = pixelization.pixel(pixelId)

        sourceCoords = SkyCoord(sourceIds["ra"] * u.degree, sourceIds["dec"] * u.degree).cartesian.xyz

        inPixel = pixelRegion.contains(*sourceCoords.value)
        if not inPixel.any():
            return sources[:0]
        pixelIds = sourceIds[inPixel]["isolated_star_id"]

        sub1 = np.clip(np.searchsorted(pixelIds, sources["isolated_star_id"]), 0, len(pixelIds) - 1)
        matched = pixelIds[sub1] == sources["isolated_star_id"]

        return sources[matched]

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        # Load specified columns from source catalogs
        names = self.collectInputNames()
        names |= {"sourceId", "coord_ra", "coord_dec"}
        for item in ["obj_index", "isolated_star_id"]:
            if item in names:
                names.remove(item)

        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(
            {"associatedSources": inputs["associatedSources"][0]}, dataId, connectionName="associatedSources"
        )

        # Loop over tract inputs, keeping only objects that in this healpix,
        # then stack in one big table.
        pixelId = dataId["healpix3"]
        associatedSourceRefs = {
            assocRef.dataId["tract"]: assocRef for assocRef in inputs["associatedSources"]
        }
        associatedSourceIdRefs = {
            assocRef.dataId["tract"]: assocRef for assocRef in inputs["associatedSourceIds"]
        }
        astrometricCorrectionRefs = {
            assocRef.dataId["tract"]: assocRef for assocRef in inputs["astrometricCorrectionCatalog"]
        }
        data = []
        for tract in associatedSourceRefs:
            tractAssociatedSources = self.loadData(associatedSourceRefs[tract], ["obj_index", "sourceId"])
            tractAssociatedSourceIds = self.loadData(
                associatedSourceIdRefs[tract], ["isolated_star_id", "ra", "dec"]
            )
            if self.config.applyAstrometricCorrections:
                astromCorrections = astrometricCorrectionRefs[tract].get(
                    parameters={"columns": self.config.astrometricCorrectionParameters.values()}
                )
            else:
                astromCorrections = None
            tractInput = {
                "associatedSources": tractAssociatedSources,
                "associatedSourceIds": tractAssociatedSourceIds,
                "astrometricCorrectionCatalog": astromCorrections,
                "sourceCatalogs": inputs["sourceCatalogs"],
                "visitTable": inputs["visitTable"],
            }
            tractData = self.callback(tractInput, dataId)
            if len(tractData) == 0:
                continue
            trimmedData = self.getHealpixOverlap(tractData, tractAssociatedSourceIds, pixelId)
            trimmedData["tract"] = tract
            data.append(trimmedData)
        data = vstack(data)

        if len(data["associatedSources"]) == 0:
            raise NoWorkFound(f"No associated sources in healpix {dataId.healpix3.id}")

        kwargs = {
            "data": data,
            "plotInfo": plotInfo,
            "camera": inputs["camera"],
        }
        outputs = self.run(**kwargs)
        self.putByBand(butlerQC, outputs, outputRefs)
