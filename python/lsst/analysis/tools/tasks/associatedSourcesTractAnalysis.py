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
from astropy.table import join, vstack
from lsst.daf.butler import DatasetProvenance
from lsst.drp.tasks.gbdesAstrometricFit import calculate_apparent_motion
from lsst.geom import Box2D
from lsst.pipe.base import NoWorkFound
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

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
            inputs["associatedSourceIds"],
            inputs["astrometricCorrectionCatalog"],
            inputs["visitTable"],
        )

    def prepareAssociatedSources(
        self,
        skymap,
        tract,
        sourceCatalogs,
        associatedSources,
        associatedSourceIds,
        astrometricCorrectionCatalog=None,
        visitTable=None,
    ):
        """Concatenate source catalogs and join on associated source IDs."""

        # Strip any provenance from tables before merging to prevent
        # warnings from conflicts being issued by astropy.utils.merge.
        for srcCat in sourceCatalogs:
            DatasetProvenance.strip_provenance_from_flat_dict(srcCat.meta)
        DatasetProvenance.strip_provenance_from_flat_dict(associatedSources.meta)
        DatasetProvenance.strip_provenance_from_flat_dict(associatedSourceIds.meta)

        # associatedSource["obj_index"] refers to the corresponding index (row)
        # in associatedSourceIds.
        index = associatedSources["obj_index"]
        associatedSources["isolated_star_id"] = associatedSourceIds["isolated_star_id"][index]

        # Keep only sources with associations
        sourceCatalogStack = vstack(sourceCatalogs, join_type="exact")
        dataJoined = join(sourceCatalogStack, associatedSources, keys="sourceId", join_type="inner")

        if astrometricCorrectionCatalog is not None:
            self.applyAstrometricCorrections(dataJoined, astrometricCorrectionCatalog, visitTable)

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

        dataWithPM = join(
            dataJoined,
            astrometricCorrectionCatalog,
            keys="isolated_star_id",
            join_type="left",
            keep_order=True,
        )

        mjds = visitTable.loc[dataWithPM["visit"]]["expMidptMJD"]
        times = astropy.time.Time(mjds, format="mjd", scale="tai")
        dataWithPM["MJD"] = times
        medianMJD = astropy.time.Time(np.median(mjds), format="mjd", scale="tai")

        raCorrection, decCorrection = calculate_apparent_motion(dataWithPM, medianMJD)

        dataJoined["coord_ra"] = dataWithPM["coord_ra"] - raCorrection.value
        dataJoined["coord_dec"] = dataWithPM["coord_dec"] - decCorrection.value

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        # Load specified columns from source catalogs
        names = self.collectInputNames()
        names |= {"sourceId", "coord_ra", "coord_dec"}
        for item in ["obj_index", "isolated_star_id"]:
            if item in names:
                names.remove(item)

        sourceCatalogs = []
        for handle in inputs["sourceCatalogs"]:
            sourceCatalogs.append(self.loadData(handle, names))
        inputs["sourceCatalogs"] = sourceCatalogs

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
