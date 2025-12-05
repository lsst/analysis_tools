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
    "SourceObjectTableAnalysisConfig",
    "SourceObjectTableAnalysisTask",
    "ObjectEpochTableConfig",
    "ObjectEpochTableTask",
)

import astropy.time
import astropy.units as u
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import numpy as np
import pandas as pd
from astropy.table import Table, join, vstack
from lsst.drp.tasks.gbdesAstrometricFit import calculate_apparent_motion
from lsst.pipe.base import AlgorithmError
from lsst.pipe.base import connectionTypes as ct
from smatch import Matcher

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class IndexMismatchError(AlgorithmError):
    """Raised if the indices in input associatedSources do not match the input
    data."""

    def __init__(self) -> None:
        super().__init__(
            "Not all sourceIds in the associated sources catalog are available in the input data."
        )

    @property
    def metadata(self) -> dict:
        """There is no metadata associated with this error."""
        return {}


class ObjectEpochTableConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "skymap"),
):
    objectCat = ct.Input(
        doc="Catalog of positions in each patch.",
        name="objectTable",
        storageClass="ArrowAstropy",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
        deferGraphConstraint=True,
    )

    epochMap = ct.Input(
        doc="Healsparse map of mean epoch of objectCat in each band.",
        name="deepCoadd_epoch_map_mean",
        storageClass="HealSparseMap",
        dimensions=("skymap", "tract", "band"),
        multiple=True,
        deferLoad=True,
    )

    objectEpochs = ct.Output(
        doc="Catalog of epochs for objectCat objects.",
        name="object_epoch",
        storageClass="ArrowAstropy",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
    )


class ObjectEpochTableConfig(pipeBase.PipelineTaskConfig, pipelineConnections=ObjectEpochTableConnections):
    bands = pexConfig.ListField(
        doc=("Bands in objectCat to be combined with `objectCat_selectors` to build objectCat column names."),
        dtype=str,
        default=["u", "g", "r", "i", "z", "y"],
    )


class ObjectEpochTableTask(pipeBase.PipelineTask):
    """Collect mean epochs for the observations that went into each object.

    TODO: DM-46202, Remove this task once the object epochs are available
    elsewhere.
    """

    ConfigClass = ObjectEpochTableConfig
    _DefaultName = "objectEpochTable"

    def getEpochs(self, cat, epochMapDict):
        """Get mean epoch of the visits corresponding to object position.

        Parameters
        ----------
        cat : `astropy.table.Table`
            Catalog containing object positions.
        epochMapDict: `dict` [`DeferredDatasetHandle`]
            Dictionary of handles for healsparse maps containing the mean epoch
            for positions in the reference catalog.

        Returns
        -------
        epochDf = `astropy.table.Table`
            Catalog with mean epoch of visits at each object position.
        """
        allEpochs = {}
        for band in self.config.bands:
            epochs = np.ones(len(cat)) * np.nan
            validPositions = np.isfinite(cat[f"{band}_ra"]) & np.isfinite(cat[f"{band}_dec"])
            if validPositions.any():
                bandEpochs = epochMapDict[band].get_values_pos(
                    cat[f"{band}_ra"][validPositions], cat[f"{band}_dec"][validPositions]
                )
                epochsValid = epochMapDict[band].get_values_pos(
                    cat[f"{band}_ra"][validPositions], cat[f"{band}_dec"][validPositions], valid_mask=True
                )
                bandEpochs[~epochsValid] = np.nan
                epochs[validPositions] = bandEpochs
            allEpochs[f"{band}_epoch"] = epochs
        allEpochs["objectId"] = cat["objectId"]

        epochTable = Table(allEpochs)
        return epochTable

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        columns = [f"{band}_{coord}" for band in self.config.bands for coord in ["ra", "dec"]]
        columns.append("objectId")

        inputs["epochMap"] = {ref.dataId["band"]: ref.get() for ref in inputs["epochMap"]}

        outputEpochRefs = {outputRef.dataId["patch"]: outputRef for outputRef in outputRefs.objectEpochs}
        for objectCatRef in inputs["objectCat"]:
            patch = objectCatRef.dataId["patch"]
            objectCat = objectCatRef.get(parameters={"columns": columns})
            epochs = self.getEpochs(objectCat, inputs["epochMap"])
            butlerQC.put(epochs, outputEpochRefs[patch])


class SourceObjectTableAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit",),
    defaultTemplates={
        "inputName": "sourceTable_visit",
        "inputCoaddName": "deep",
        "associatedSourcesInputName": "isolated_star_presources",
        "associatedSourceIdsInputName": "isolated_star_presource_associations",
        "outputName": "sourceObjectTable",
    },
):
    data = ct.Input(
        doc="Visit based source table to load from the butler",
        name="sourceTable_visit",
        storageClass="ArrowAstropy",
        dimensions=("visit",),
        deferLoad=True,
    )

    associatedSources = ct.Input(
        doc="Table of associated sources",
        name="{associatedSourcesInputName}",
        storageClass="ArrowAstropy",
        multiple=True,
        deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
        deferGraphConstraint=True,
    )

    associatedSourceIds = ct.Input(
        doc="Table containing unique ids for the associated sources",
        name="{associatedSourceIdsInputName}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        multiple=True,
        dimensions=("instrument", "skymap", "tract"),
        deferGraphConstraint=True,
    )

    refCat = ct.Input(
        doc="Catalog of positions to use as reference.",
        name="objectTable",
        storageClass="DataFrame",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
        deferGraphConstraint=True,
    )
    astrometricCorrectionCatalog = ct.Input(
        doc="Catalog containing proper motions and parallaxes.",
        name="isolated_star_stellar_motions",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "skymap", "tract"),
        multiple=True,
        deferLoad=True,
    )
    refCatEpochs = ct.Input(
        doc="Catalog of epochs for refCat objects.",
        name="object_epoch",
        storageClass="ArrowAstropy",
        dimensions=["skymap", "tract", "patch"],
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
            self.inputs.remove("astrometricCorrectionCatalog")
            self.inputs.remove("refCatEpochs")
            self.inputs.remove("visitTable")


class SourceObjectTableAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=SourceObjectTableAnalysisConnections
):
    ra_column = pexConfig.Field(
        doc="Name of column in refCat to use for right ascension.",
        dtype=str,
        default="r_ra",
    )
    dec_column = pexConfig.Field(
        doc="Name of column in refCat to use for declination.",
        dtype=str,
        default="r_dec",
    )
    epoch_column = pexConfig.Field(
        doc=(
            "Name of column in refCat corresponding to the epoch to which "
            "sources will be shifted. Should correspond to the positions in "
            "`ra_column` and `dec_column`."
        ),
        dtype=str,
        default="r_epoch",
    )
    refCat_bands = pexConfig.ListField(
        doc=("Bands in refCat to be combined with `refCat_selectors` to build refCat column names."),
        dtype=str,
        default=["u", "g", "r", "i", "z", "y"],
    )
    refCat_selectors = pexConfig.ListField(
        doc=(
            "Remove objects for which these flags are true. These strings are combined with `refCat_bands`"
            " to build the full refCat column names"
        ),
        dtype=str,
        default=["pixelFlags_saturated", "pixelFlags_saturatedCenter"],
    )
    refCatMatchingRadius = pexConfig.Field(
        dtype=float,
        default=1.0,
        doc=(
            "Radius in mas with which to match the mean positions of the sources with the positions in the"
            " reference catalog."
        ),
    )
    applyAstrometricCorrections = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Apply proper motions and parallaxes to source positions.",
    )
    correctionsMatchingRadius = pexConfig.Field(
        dtype=float,
        default=0.2,
        doc=(
            "Radius in mas with which to match the mean positions of the sources with the positions in the"
            " astrometricCorrectionCatalog."
        ),
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

    def setDefaults(self):
        super().setDefaults()
        from ..atools import TargetRefCatDeltaColorMetrics

        self.atools.astromColorDiffMetrics = TargetRefCatDeltaColorMetrics


class SourceObjectTableAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SourceObjectTableAnalysisConfig
    _DefaultName = "sourceObjectTableAnalysis"

    def callback(self, inputs, dataId):
        """Callback function to be used with reconstructor."""
        return self.prepareAssociatedSources(
            dataId["visit"],
            inputs["data"],
            inputs["associatedSources"],
            inputs["associatedSourceIds"],
            inputs["refCat"],
            inputs["visitTable"],
            inputs["astrometricCorrectionCatalog"],
        )

    def applyAstrometricCorrections(
        self, isolatedSources, astrometricCorrectionCatalog, visitTable, visit, refEpochs
    ):
        """Shift source positions to match the epoch of the reference catalog
        objects.

        Parameters
        ----------
        isolatedSources : `astropy.table.Table`
            Catalog of sources which will be modified in place with the
            astrometric corrections.
        astrometricCorrectionCatalog : `astropy.table.Table`
            Catalog with proper motion and parallax information.
        visitTable : `pd.DataFrame`
            Catalog containing the epoch for the visit corresponding to the
            isolatedSources.
        visit : `int`
            Identifier of the isolatedSources' visit.
        """
        if visitTable.index.name is None:
            # The expected index may or may not be set, depending on whether
            # the table was written originally as a DataFrame or something else
            # Parquet-friendly.
            visitTable.set_index("visitId", inplace=True)
        sourceMjd = visitTable.loc[visit]["expMidptMJD"]

        # Get target date from reference catalog
        targetEpochs = refEpochs.to_numpy()
        # There may not be a valid reference epoch on the edge of a given
        # region. Do not make an astrometric correction for any sources on the
        # edge.
        targetEpochs[~np.isfinite(targetEpochs)] = sourceMjd

        # Get the stellar motion catalog into the right format:
        for key, value in self.config.astrometricCorrectionParameters.items():
            astrometricCorrectionCatalog.rename_column(value, key)
        astrometricCorrectionCatalog["ra"] *= u.degree
        astrometricCorrectionCatalog["dec"] *= u.degree
        astrometricCorrectionCatalog["pmRA"] *= u.mas / u.yr
        astrometricCorrectionCatalog["pmDec"] *= u.mas / u.yr
        astrometricCorrectionCatalog["parallax"] *= u.mas

        joinedData = join(
            isolatedSources[["isolated_star_id"]],
            astrometricCorrectionCatalog,
            keys="isolated_star_id",
            join_type="left",
            keep_order=True,
            metadata_conflicts="silent",
        )
        joinedData["MJD"] = astropy.time.Time(sourceMjd, format="mjd", scale="tai")

        raCorrection, decCorrection = calculate_apparent_motion(
            joinedData, astropy.time.Time(targetEpochs, format="mjd", scale="tai")
        )

        isolatedSources["coord_ra"] -= raCorrection.value
        isolatedSources["coord_dec"] -= decCorrection.value

    def prepareAssociatedSources(
        self,
        visit,
        data,
        associatedSourceRefs,
        associatedSourceIdRefs,
        refCats,
        visitTable,
        astrometricCorrectionCatalog,
    ):
        """Match isolated sources with reference objects and shift the sources
        to the object epochs if `self.config.applyAstrometricCorrections` is
        True.

        Parameters
        ----------
        visit : `int`
            Identifier of the visit corresponding to the data.
        data : `astropy.table.Table`
            Catalog of sources to be associated.
        associatedSourceRefs : `list` [`DeferredDatasetHandle`]
            Handle for the catalogs of isolated sources. There will be multiple
            if the visit overlaps with multiple tracts.
        refCats : `list` [`pd.DataFrame`]
            Catalog of objects with which the sources will be compared.
        visitTable : `pd.DataFrame`
            Catalog containing the epoch for the visit corresponding to the
            isolatedSources.
        astrometricCorrectionCatalog : `astropy.table.Table`
            Catalog with proper motion and parallax information.
        """
        isolatedSources = []
        associatedSourceIds = {
            ref.dataId["tract"]: ref.get(parameters={"columns": ["isolated_star_id"]})
            for ref in associatedSourceIdRefs
        }
        for associatedSourceRef in associatedSourceRefs:
            tract = associatedSourceRef.dataId["tract"]
            associatedSources = associatedSourceRef.get(
                parameters={"columns": ["visit", "sourceId", "obj_index"]}
            )
            index = associatedSources["obj_index"]
            associatedSources["isolated_star_id"] = associatedSourceIds[tract]["isolated_star_id"][index]

            visit_sources = associatedSources[associatedSources["visit"] == visit]
            try:
                visitData = data.loc[visit_sources["sourceId"]]
                visitData["isolated_star_id"] = visit_sources["isolated_star_id"]
                isolatedSources.append(visitData)
            except KeyError:
                raise IndexMismatchError()
        isolatedSources = vstack(isolatedSources)

        if len(isolatedSources) == 0:
            raise pipeBase.NoWorkFound(f"No isolated sources found for visit {visit}")

        with Matcher(np.asarray(isolatedSources["coord_ra"]), np.asarray(isolatedSources["coord_dec"])) as m:
            idx, isolatedMatchIndices, refMatchIndices, dists = m.query_radius(
                np.asarray(refCats[self.config.ra_column]),
                np.asarray(refCats[self.config.dec_column]),
                self.config.refCatMatchingRadius / 3600.0,
                return_indices=True,
            )

        matchIS = isolatedSources[isolatedMatchIndices]
        # Apply proper motions and parallaxes to visit sources.
        if self.config.applyAstrometricCorrections:
            refCatEpochs = refCats[self.config.epoch_column].iloc[refMatchIndices]
            self.applyAstrometricCorrections(
                matchIS, astrometricCorrectionCatalog, visitTable, visit, refCatEpochs
            )

        matchRef = refCats.iloc[refMatchIndices]
        matchIS = matchIS.to_pandas()

        allCat = pd.concat([matchRef.reset_index(), matchIS.reset_index()], axis=1)
        return allCat

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)

        # Get isolated sources:
        visit = inputs["data"].dataId["visit"]
        band = inputs["data"].dataId["band"]
        names = self.collectInputNames()
        names -= {self.config.ra_column, self.config.dec_column}
        names.add("sourceId")
        data = inputs["data"].get(parameters={"columns": names})
        data.add_index("sourceId")
        inputs["data"] = data

        if self.config.applyAstrometricCorrections:
            refCatEpochs = {
                epochTable.dataId["patch"]: epochTable.get() for epochTable in inputs["refCatEpochs"]
            }
        # Get objects:
        allRefCats = []
        refCatSelectors = [
            f"{refCatBand}_{selector}"
            for refCatBand in self.config.refCat_bands
            for selector in self.config.refCat_selectors
        ]

        for refCatRef in inputs["refCat"]:
            refCat = refCatRef.get(
                parameters={
                    "columns": ["detect_isPrimary", self.config.ra_column, self.config.dec_column, "objectId"]
                    + refCatSelectors
                }
            )
            refCat.set_index("objectId")
            if self.config.applyAstrometricCorrections:
                refCat = pd.merge(refCat, refCatEpochs[refCatRef.dataId["patch"]].to_pandas(), on="objectId")
            goodInds = (
                refCat["detect_isPrimary"]
                & np.isfinite(refCat[self.config.ra_column])
                & np.isfinite(refCat[self.config.dec_column])
            )
            goodInds &= ~refCat[refCatSelectors].any(axis=1)
            allRefCats.append(refCat[goodInds])

        refCat = pd.concat(allRefCats)
        inputs["refCat"] = refCat
        if len(refCat) == 0:
            raise pipeBase.NoWorkFound(f"No reference catalog objects found to associate with visit {visit}")

        if self.config.applyAstrometricCorrections:
            pmCats = []
            for astrometricCorrectionCatalogRef in inputs["astrometricCorrectionCatalog"]:
                pmCat = astrometricCorrectionCatalogRef.get(
                    parameters={"columns": self.config.astrometricCorrectionParameters.values()}
                )
                pmCats.append(pmCat)
            inputs["astrometricCorrectionCatalog"] = vstack(pmCats, metadata_conflicts="silent")
        else:
            inputs["astrometricCorrectionCatalog"] = None
            inputs["visitTable"] = None

        allCat = self.callback(inputs, dataId)

        outputs = self.run(data=allCat, bands=band, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)
