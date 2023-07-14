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

__all__ = ("CatalogMatchConfig", "CatalogMatchTask")


import numpy as np
from astropy.time import Time
from astropy.table import Table, hstack
from smatch import Matcher
import os.path

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask
from lsst.skymap import BaseSkyMap
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
import lsst.geom
from ..actions.vector import (
    CoaddPlotFlagSelector,
    GalaxySelector,
    SnSelector,
    StarSelector,
    VisitPlotFlagSelector,
)
from ..interfaces import VectorAction


class CatalogMatchConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"targetCatalog": "objectTable_tract", "refCatalog": "ps1_pv3_3pi_20170110"},
):
    catalog = pipeBase.connectionTypes.Input(
        doc="The tract-wide catalog to make plots from.",
        storageClass="ArrowAstropy",
        name="{targetCatalog}",
        dimensions=("tract", "skymap"),
        deferLoad=True,
    )

    refCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The reference catalog to match to loaded input catalog sources.",
        name="{refCatalog}",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    skymap = pipeBase.connectionTypes.Input(
        doc="The skymap for the tract",
        storageClass="SkyMap",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap",),
    )

    matchedCatalog = pipeBase.connectionTypes.Output(
        doc="Catalog with matched target and reference objects with separations",
        name="{targetCatalog}_{refCatalog}_match",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )


class CatalogMatchConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=CatalogMatchConnections
):

    referenceCatalogLoader = pexConfig.ConfigurableField(
        target=LoadReferenceCatalogTask,
        doc="Reference catalog loader",
    )

    epoch = pexConfig.Field[float](
        doc="Epoch to which reference objects are shifted.",
        default=2015.0,
    )

    filterNames = pexConfig.ListField[str](
        doc="Physical filter names to persist downstream.",
        default=["u", "g", "r", "i", "z", "y"],
    )

    selectorBands = pexConfig.ListField[str](
        doc="Band to use when selecting objects, primarily for extendedness.",
        default=["i"],
    )

    selectorActions = ConfigurableActionStructField[VectorAction](
        doc="Which selectors to use to narrow down the data for QA plotting.",
        default={"flagSelector": CoaddPlotFlagSelector()},
    )

    sourceSelectorActions = ConfigurableActionStructField[VectorAction](
        doc="What types of sources to use.",
        default={"sourceSelector": StarSelector()},
    )

    extraColumnSelectors = ConfigurableActionStructField[VectorAction](
        doc="Other selectors that are not used in this task, but whose columns" "may be needed downstream",
        default={"selector1": SnSelector(), "selector2": GalaxySelector()},
    )

    extraColumns = pexConfig.ListField[str](
        doc="Other catalog columns to persist to downstream tasks",
        default=["x", "y", "patch", "ebv"],
    )

    extraPerBandColumns = pexConfig.ListField[str](
        doc="Other columns to load that should be loaded for each band individually.",
        default=["cModelFlux"],
    )

    matchRadius = pexConfig.Field[float](
        doc="The radius to use for matching, in arcsecs.",
        default=1.0,
    )

    raColumn = pexConfig.Field[str](doc="RA column.", default="coord_ra")
    decColumn = pexConfig.Field[str](doc="Dec column.", default="coord_dec")
    patchColumn = pexConfig.Field[str](doc="Patch column.", default="patch")

    # in the validation, we need to confirm that all of the filter names are in the filter map.

    # do a setDefaults and set however you like.
    def setDefaults(self):
        super().setDefaults()
        self.referenceCatalogLoader.doReferenceSelection = False
        self.referenceCatalogLoader.doApplyColorTerms = False

class CatalogMatchTask(pipeBase.PipelineTask):
    """The base task for matching catalogs. Figures out which columns
    it needs to grab for the downstream tasks and then matches the
    two tables together and returns the matched and joined table
    including the extra columns.
    """
    ConfigClass = CatalogMatchConfig
    _DefaultName = "analysisToolsCatalogMatch"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Implemented in the inherited tasks"""
        pass

    def run(self, *, catalog, loadedRefCat, bands):
        """Takes the two catalogs and returns the matched one.

        Parameters
        ----------
        `catalog` : astropy.table.Table
            The catalog to be matched
        `loadedRefCat` : astropy.table.Table
            The loaded reference catalog
        `bands` : list
            A list of bands to apply the selectors in

        Returns
        -------
        `matchedCatalog` : astropy.table.Table

        Notes
        -----
        Performs an RA/Dec match that returns the closest match
        within the match radius which defaults to 1.0 arcsecond.
        Applies the suffix, _target, to the catalog being matched
        and _ref to the reference catalog.
        """
        # Apply the selectors to the catalog
        mask = np.ones(len(catalog), dtype=bool)
        for selector in self.config.sourceSelectorActions:
            for band in self.config.selectorBands:
                mask &= selector(catalog, band=band).astype(bool)

        targetCatalog = catalog[mask]

        if (len(targetCatalog) == 0) or (len(loadedRefCat) == 0):
            refMatchIndices = np.array([], dtype=np.int64)
            targetMatchIndices = np.array([], dtype=np.int64)
            dists = np.array([], dtype=np.float64)
        else:
            # Run the matcher.

            # This all assumes that everything is in degrees.
            # Which I think is okay, but the current task allows different units.
            # Need to configure match radius, either in this task or a subtask.
            with Matcher(loadedRefCat["ra"], loadedRefCat["dec"]) as m:
                idx, refMatchIndices, targetMatchIndices, dists = m.query_radius(
                    targetCatalog[self.config.raColumn],
                    targetCatalog[self.config.decColumn],
                    self.config.matchRadius / 3600.0,
                    return_indices=True,
                )

            # Convert degrees to arcseconds.
            dists *= 3600.0

        targetCatalogMatched = targetCatalog[targetMatchIndices]
        loadedRefCatMatched = loadedRefCat[refMatchIndices]

        targetCols = targetCatalogMatched.columns.copy()
        for col in targetCols:
            targetCatalogMatched.rename_column(col, col + "_target")
        refCols = loadedRefCatMatched.columns.copy()
        for col in refCols:
            loadedRefCatMatched.rename_column(col, col + "_ref")

        for (i, band) in enumerate(bands):
            loadedRefCatMatched[band + "_mag_ref"] = loadedRefCatMatched["refMag_ref"][:, i]
            loadedRefCatMatched[band + "_magErr_ref"] = loadedRefCatMatched["refMagErr_ref"][:, i]
        loadedRefCatMatched.remove_column("refMag_ref")
        loadedRefCatMatched.remove_column("refMagErr_ref")
        tMatched = hstack([targetCatalogMatched, loadedRefCatMatched], join_type="exact")
        tMatched["matchDistance"] = dists

        return pipeBase.Struct(matchedCatalog=tMatched)

    def prepColumns(self, bands):
        """Get all the columns needed for downstream tasks.
        Both those from the selectors and those specified in the
        config options.
        """

        bandColumns = []
        for band in bands:
            for col in self.config.extraPerBandColumns:
                bandColumns.append(band + "_" + col)

        columns = [
            self.config.raColumn,
            self.config.decColumn,
        ] + self.config.extraColumns.list() + bandColumns

        if self.config.patchColumn is not "":
            columns.append(self.config.patchColumn)

        selectorBands = list(set(list(bands) + self.config.selectorBands.list()))
        for selectorAction in [
            self.config.selectorActions,
            self.config.sourceSelectorActions,
            self.config.extraColumnSelectors,
        ]:
            for selector in selectorAction:
                for band in selectorBands:
                    selectorSchema = selector.getFormattedInputSchema(band=band)
                    columns += [s[0] for s in selectorSchema]

        return columns

    def _loadRefCat(self, loaderTask, tractInfo):
        """Load the reference catalog that covers the
        catalog that is to be matched to.

        Parameters
        ----------
        `loaderTask` : lsst.pipe.tasks.loadReferenceCatalog.loadReferenceCatalogTask
        `tractInfo` : lsst.skymap.tractInfo.ExplicitTractInfo
            The tract information to get the sky location from

        Returns
        -------
        `loadedRefCat` : astropy.table.Table
            The reference catalog that covers the input catalog.
        """
        boundingCircle = tractInfo.getOuterSkyPolygon().getBoundingCircle()
        center = lsst.geom.SpherePoint(boundingCircle.getCenter())
        radius = boundingCircle.getOpeningAngle()

        epoch = Time(self.config.epoch, format="decimalyear")

        # This is always going to return degrees.
        loadedRefCat = loaderTask.getSkyCircleCatalog(center, radius, self.config.filterNames, epoch=epoch)

        return Table(loadedRefCat)
