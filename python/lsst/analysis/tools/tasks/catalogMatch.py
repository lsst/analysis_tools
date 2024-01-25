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


import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import numpy as np
from astropy.table import Table, hstack, vstack
from astropy.time import Time
from lsst.pex.config.configurableActions import ConfigurableActionStructField
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask
from lsst.skymap import BaseSkyMap
from smatch import Matcher

from ..actions.vector import (
    CoaddPlotFlagSelector,
    GalaxySelector,
    MatchingFlagSelector,
    SnSelector,
    StarSelector,
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


class CatalogMatchConfig(pipeBase.PipelineTaskConfig, pipelineConnections=CatalogMatchConnections):
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
        default={"flagSelector": MatchingFlagSelector()},
    )

    sourceSelectorActions = ConfigurableActionStructField[VectorAction](
        doc="What types of sources to use.",
        default={},
    )

    extraColumnSelectors = ConfigurableActionStructField[VectorAction](
        doc="Other selectors that are not used in this task, but whose columns" "may be needed downstream",
        default={
            "selector1": SnSelector(),
            "selector2": StarSelector(),
            "selector3": GalaxySelector(),
            "selector4": CoaddPlotFlagSelector(),
        },
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

    targetRaColumn = pexConfig.Field[str](
        doc="RA column name for the target (being matched) catalog.",
        default="coord_ra",
    )

    targetDecColumn = pexConfig.Field[str](
        doc="Dec column name for the target (being matched) catalog.",
        default="coord_dec",
    )

    refRaColumn = pexConfig.Field[str](
        doc="RA column name for the reference (being matched to) catalog.",
        default="ra",
    )

    refDecColumn = pexConfig.Field[str](
        doc="Dec column name for the reference (being matched to) catalog.",
        default="dec",
    )

    raColumn = pexConfig.Field[str](
        doc="RA column.",
        default="coord_ra",
        deprecated="This field was replaced with targetRaColumn and is unused. Will be removed after v27.",
    )

    decColumn = pexConfig.Field[str](
        doc="Dec column.",
        default="coord_dec",
        deprecated="This field was replaced with targetDecColumn and is unused. Will be removed after v27.",
    )

    patchColumn = pexConfig.Field[str](doc="Patch column.", default="patch")

    matchesRefCat = pexConfig.Field[bool](
        doc="Is the catalog being matched to stored as a reference catalog?",
        default=False,
    )

    returnNonMatches = pexConfig.Field[bool](
        doc="Return the rows of the reference catalog that didn't get matched?",
        default=False,
    )

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

    def run(self, *, targetCatalog, refCatalog, bands):
        """Takes the two catalogs and returns the matched one.

        Parameters
        ----------
        `targetCatalog` : astropy.table.Table
            The catalog to be matched
        `refCatalog` : astropy.table.Table
            The catalog to be matched to
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
        and _ref to the reference catalog being matched to.
        """
        # Apply the selectors to the catalog
        mask = np.ones(len(targetCatalog), dtype=bool)
        for selector in self.config.sourceSelectorActions:
            for band in self.config.selectorBands:
                mask &= selector(targetCatalog, band=band).astype(bool)

        targetCatalog = targetCatalog[mask]

        if (len(targetCatalog) == 0) or (len(refCatalog)) == 0:
            refMatchIndices = np.array([], dtype=np.int64)
            targetMatchIndices = np.array([], dtype=np.int64)
            dists = np.array([], dtype=np.float64)
        else:
            # Run the matcher.

            # This all assumes that everything is in degrees.
            # Which I think is okay, but the current task
            # allows different units. Need to configure match
            # radius, either in this task or a subtask.

            # Get rid of entries in the refCat with non-finite RA/Dec values.
            refRas = refCatalog[self.config.refRaColumn]
            refDecs = refCatalog[self.config.refDecColumn]
            refRaDecFiniteMask = np.isfinite(refRas) & np.isfinite(refDecs)
            refCatalog = refCatalog[refRaDecFiniteMask]
            with Matcher(refCatalog[self.config.refRaColumn], refCatalog[self.config.refDecColumn]) as m:
                idx, refMatchIndices, targetMatchIndices, dists = m.query_radius(
                    targetCatalog[self.config.targetRaColumn],
                    targetCatalog[self.config.targetDecColumn],
                    self.config.matchRadius / 3600.0,
                    return_indices=True,
                )

            # Convert degrees to arcseconds.
            dists *= 3600.0

        targetCatalogMatched = targetCatalog[targetMatchIndices]
        refCatalogMatched = refCatalog[refMatchIndices]

        targetCols = targetCatalogMatched.columns.copy()
        for col in targetCols:
            targetCatalogMatched.rename_column(col, col + "_target")
        refCols = refCatalogMatched.columns.copy()
        for col in refCols:
            refCatalogMatched.rename_column(col, col + "_ref")

        if self.config.returnNonMatches:
            unmatchedIndices = list(set(np.arange(0, len(refCatalog))) - set(refMatchIndices))
            refCatalogNotMatched = refCatalog[unmatchedIndices]
            # We need to set the relevant flag columns to
            # true or false so that they make it through the
            # selectors even though the none matched sources
            # don't have values for those columns.
            trueFlagCols = []
            falseFlagCols = []
            for selectorAction in [self.config.selectorActions, self.config.extraColumnSelectors]:
                for selector in selectorAction:
                    try:
                        for flag in selector.selectWhenTrue:
                            trueFlagCols.append(flag)
                        for flag in selector.selectWhenFalse:
                            falseFlagCols.append(flag)
                    except AttributeError:
                        continue
            for col in refCols:
                refCatalogNotMatched.rename_column(col, col + "_ref")
            for col in targetCols:
                refCatalogNotMatched[col] = [np.nan] * len(refCatalogNotMatched)
            for col in trueFlagCols:
                refCatalogNotMatched[col] = [True] * len(refCatalogNotMatched)
            for col in falseFlagCols:
                refCatalogNotMatched[col] = [False] * len(refCatalogNotMatched)

        if self.config.matchesRefCat:
            for i, band in enumerate(bands):
                refCatalogMatched[band + "_mag_ref"] = refCatalogMatched["refMag_ref"][:, i]
                refCatalogMatched[band + "_magErr_ref"] = refCatalogMatched["refMagErr_ref"][:, i]
            refCatalogMatched.remove_column("refMag_ref")
            refCatalogMatched.remove_column("refMagErr_ref")

            if self.config.returnNonMatches:
                for i, band in enumerate(bands):
                    refCatalogNotMatched[band + "_mag_ref"] = refCatalogNotMatched["refMag_ref"][:, i]
                    refCatalogNotMatched[band + "_magErr_ref"] = refCatalogNotMatched["refMagErr_ref"][:, i]
                refCatalogNotMatched.remove_column("refMag_ref")
                refCatalogNotMatched.remove_column("refMagErr_ref")

        tMatched = hstack([targetCatalogMatched, refCatalogMatched], join_type="exact")
        tMatched["matchDistance"] = dists

        if self.config.returnNonMatches:
            refCatalogNotMatched["matchDistance"] = [np.nan] * len(refCatalogNotMatched)
            tMatched = vstack([tMatched, refCatalogNotMatched])

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

        columns = (
            [
                self.config.targetRaColumn,
                self.config.targetDecColumn,
            ]
            + self.config.extraColumns.list()
            + bandColumns
        )

        if self.config.patchColumn != "":
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
        `loaderTask` :
            lsst.pipe.tasks.loadReferenceCatalog.loadReferenceCatalogTask
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
        try:
            loadedRefCat = loaderTask.getSkyCircleCatalog(
                center, radius, self.config.filterNames, epoch=epoch
            )
        except RuntimeError as e:
            raise pipeBase.NoWorkFound(e)

        return Table(loadedRefCat)
