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

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import numpy as np
import pandas as pd
from astropy.table import vstack
from lsst.drp.tasks.gbdesAstrometricFit import calculate_apparent_motion
from lsst.pipe.base import connectionTypes as ct
from smatch import Matcher

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class ObjectEpochTableConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={
        # "inputName": "sourceTable_visit",
        # "inputCoaddName": "deep",
        # "associatedSourcesInputName": "isolated_star_sources",
        # "outputName": "sourceObjectTable",
    },
):
    refCat = ct.Input(
        doc="Catalog of positions to use as reference.",
        name="objectTable",
        storageClass="DataFrame",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
        deferGraphConstraint=True,
    )

    refCatEpochMap = ct.Input(
        doc="Healsparse map of mean epoch of refCat.",
        name="deepCoadd_epoch_map_mean",
        storageClass="HealSparseMap",
        dimensions=("skymap", "tract", "band"),
        multiple=True,
        deferLoad=True,
    )

    objectEpochs = ct.Output(
        doc="Catalog of epochs for refCat objects.",
        name="objectEpochTable",
        storageClass="DataFrame",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
    )


class ObjectEpochTableConfig(pipeBase.PipelineTaskConfig, pipelineConnections=ObjectEpochTableConnections):
    refCatBands = pexConfig.ListField(
        doc=("Bands in refCat to be combined with `refCat_selectors` to build refCat column names."),
        dtype=str,
        default=["u", "g", "r", "i", "z", "y"],
    )


class ObjectEpochTableTask(pipeBase.PipelineTask):
    """Collect mean epochs for the observations that went into each object."""

    ConfigClass = ObjectEpochTableConfig
    _DefaultName = "objectEpochTable"

    def getEpochs(self, cat, epochMapDict):
        """Get mean epoch of the visits corresponding to object position.

        Parameters
        ----------
        cat : `pd.DataFrame`
            Catalog containing object positions.
        refCatEpochMap: `dict` [`DeferredDatasetHandle`]
            Dictionary of handles for healsparse maps containing the mean epoch
            for positions in the reference catalog.

        Returns
        -------
        epochDf = `pd.DataFrame`
            Catalog with mean epoch of visits at each object positions.
        """
        epochs = {}
        for band in self.config.refCatBands:
            refEpochs = np.ones(len(cat)) * np.nan
            validPositions = np.isfinite(cat[f"{band}_ra"]) & np.isfinite(cat[f"{band}_dec"])
            if validPositions.any():
                bandEpochs = epochMapDict[band].get_values_pos(
                    cat[f"{band}_ra"][validPositions], cat[f"{band}_dec"][validPositions]
                )
                refEpochsValid = epochMapDict[band].get_values_pos(
                    cat[f"{band}_ra"][validPositions], cat[f"{band}_dec"][validPositions], valid_mask=True
                )
                bandEpochs[~refEpochsValid] = np.nan
                refEpochs[validPositions] = bandEpochs
            epochs[f"{band}_epoch"] = refEpochs
        epochs["objectId"] = cat.index

        epochDf = pd.DataFrame(epochs).set_index("objectId")
        return epochDf

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        columns = [f"{band}_{coord}" for band in self.config.refCatBands for coord in ["ra", "dec"]]

        inputs["refCatEpochMap"] = {ref.dataId["band"]: ref.get() for ref in inputs["refCatEpochMap"]}

        outputEpochRefs = {outputRef.dataId["patch"]: outputRef for outputRef in outputRefs.objectEpochs}
        for refCatRef in inputs["refCat"]:
            patch = refCatRef.dataId["patch"]
            refCat = refCatRef.get(parameters={"columns": columns})
            refEpochs = self.getEpochs(refCat, inputs["refCatEpochMap"])
            butlerQC.put(refEpochs, outputEpochRefs[patch])


class SourceObjectTableAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit",),
    defaultTemplates={
        "inputName": "sourceTable_visit",
        "inputCoaddName": "deep",
        "associatedSourcesInputName": "isolated_star_sources",
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

    refCat = ct.Input(
        doc="Catalog of positions to use as reference.",
        name="objectTable",
        storageClass="DataFrame",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
        deferGraphConstraint=True,
    )
    properMotionCatalog = ct.Input(
        doc="Catalog containing proper motions.",
        name="gbdesAstrometricFit_starCatalog",
        storageClass="ArrowNumpyDict",
        dimensions=("instrument", "skymap", "tract", "physical_filter"),
        multiple=True,
        deferLoad=True,
    )
    refCatEpochs = ct.Input(
        doc="Catalog of epochs for refCat objects.",
        name="objectEpochTable",
        storageClass="DataFrame",
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

        if not config.applyProperMotions:
            self.inputs.remove("properMotionCatalog")
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
        doc="Name of column in refCat to use for epoch.",
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
    applyProperMotions = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Apply proper motions to source positions.",
    )


class SourceObjectTableAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SourceObjectTableAnalysisConfig
    _DefaultName = "sourceObjectTableAnalysis"

    def callback(self, inputs, dataId):
        """Callback function to be used with reconstructor."""
        return self.prepareAssociatedSources(
            dataId["visit"],
            inputs["data"],
            inputs["associatedSources"],
            inputs["refCat"],
            inputs["visitTable"],
            inputs["properMotionCatalog"],
        )

    def applyProperMotionCorrections(
        self, isolatedSources, properMotionCatalog, visitTable, visit, refEpochs
    ):
        """Shift source positions to match the epoch of the reference catalog
        objects.

        Parameters
        ----------
        isolatedSources : `astropy.table.Table`
            Catalog of sources to which proper motion corrections will be made.
        properMotionCatalog : `pd.DataFrame`
            Catalog with proper motion and parallax information.
        visitTable : `pd.DataFrame`
            Catalog containing the epoch for the visit corresponding to the
            isolatedSources.
        visit : `int`
            Identifier of the isolatedSources' visit.
        """
        sourceMjd = visitTable.loc[visit]["expMidptMJD"]

        # Get target date from reference catalog
        targetEpochs = refEpochs.to_numpy()
        # There may not be a valid reference epoch on the edge of a given
        # region. Do not make a proper motion correction for any sources on the
        # edge.
        targetEpochs[~np.isfinite(targetEpochs)] = sourceMjd

        with Matcher(isolatedSources["coord_ra"], isolatedSources["coord_dec"]) as m:
            idx, i1, i2, d = m.query_radius(
                properMotionCatalog["starX"], properMotionCatalog["starY"], 0.2 / 3600, return_indices=True
            )
        pmDf = properMotionCatalog.rename(
            columns={
                "starX": "ra",
                "starY": "dec",
                "starPMx": "pmRA",
                "starPMy": "pmDec",
                "starParallax": "parallax",
            }
        )
        pmDf.insert(5, "MJD", sourceMjd)
        raCorrection, decCorrection = calculate_apparent_motion(pmDf.iloc[i2], targetEpochs[i1])
        isolatedSources["coord_ra"][i1] -= raCorrection.value
        isolatedSources["coord_dec"][i1] -= decCorrection.value
        return isolatedSources

    def prepareAssociatedSources(
        self, visit, data, associatedSourceRefs, refCats, visitTable, properMotionCatalog
    ):
        """Match isolated sources with reference objects and shift the sources
        to the object epochs if the properMotionCatalog is not None.

        Parameters
        ----------
        visit : `int`
            Identifier of the visit corresponding to the data.
        data : `astropy.table.Table`
            Catalog of sources.
        associatedSourceRefs : `list` [`DeferredDatasetHandle`]
            Handle for the catalogs of isolated sources. There will be multiple
            if the visit overlaps with multiple tracts.
        refCats : `list` [`pd.DataFrame`]
            Catalog of objects with which the sources will be compared.
        visitTable : `pd.DataFrame`
            Catalog containing the epoch for the visit corresponding to the
            isolatedSources.
        properMotionCatalog : `pd.DataFrame`
            Catalog with proper motion and parallax information.
        """
        isolatedSources = []
        for associatedSourceRef in associatedSourceRefs:
            associatedSources = associatedSourceRef.get(parameters={"columns": ["visit", "source_row"]})
            visit_sources = associatedSources[associatedSources["visit"] == visit]
            isolatedSources.append(data[visit_sources["source_row"]])
        isolatedSources = vstack(isolatedSources)

        if len(isolatedSources) == 0:
            raise pipeBase.NoWorkFound(f"No isolated sources found for visit {visit}")

        with Matcher(isolatedSources["coord_ra"], isolatedSources["coord_dec"]) as m:
            idx, isolatedMatchIndices, refMatchIndices, dists = m.query_radius(
                refCats[self.config.ra_column].values,
                refCats[self.config.dec_column].values,
                1 / 3600.0,
                return_indices=True,
            )

        refCatEpochs = refCats[self.config.epoch_column].iloc[refMatchIndices]
        matchIS = isolatedSources[isolatedMatchIndices]
        # Get proper motions for visit sources
        if self.config.applyProperMotions:
            matchIS = self.applyProperMotionCorrections(
                matchIS, properMotionCatalog, visitTable, visit, refCatEpochs
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
        data = inputs["data"].get(parameters={"columns": names})
        inputs["data"] = data

        if self.config.applyProperMotions:
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
                    "columns": ["detect_isPrimary", self.config.ra_column, self.config.dec_column]
                    + refCatSelectors
                }
            )
            if self.config.applyProperMotions:
                refCat = refCat.join(refCatEpochs[refCatRef.dataId["patch"]])
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

        if self.config.applyProperMotions:
            pmCats = []
            for properMotionCatalogRef in inputs["properMotionCatalog"]:
                pmCat = properMotionCatalogRef.get(
                    parameters={"columns": ["starX", "starY", "starPMx", "starPMy", "starParallax"]}
                )
                pmCats.append(pd.DataFrame(pmCat))
            inputs["properMotionCatalog"] = pd.concat(pmCats)
        else:
            inputs["properMotionCatalog"] = None
            inputs["visitTable"] = None

        allCat = self.callback(inputs, dataId)

        outputs = self.run(data=allCat, bands=band, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)
