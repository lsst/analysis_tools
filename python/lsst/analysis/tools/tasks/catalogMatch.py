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

__all__ = ("CatalogMatchConfig", "CatalogMatchTask", "AstropyMatchConfig", "AstropyMatchTask")

import astropy.units as units
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from lsst.meas.algorithms import ReferenceObjectLoader
from lsst.pex.config.configurableActions import ConfigurableActionStructField
from lsst.skymap import BaseSkyMap

from ..actions.vector import (
    CoaddPlotFlagSelector,
    GalaxySelector,
    SnSelector,
    StarSelector,
    VisitPlotFlagSelector,
)
from ..interfaces import VectorAction


class AstropyMatchConfig(pexConfig.Config):
    maxDistance = pexConfig.Field[float](
        doc="Max distance between matches in arcsec",
        default=1.0,
    )
    refCatUnits = pexConfig.Field[str](
        doc="Units of the reference catalog coordinates",
        default="degree",
    )
    targetCatUnits = pexConfig.Field[str](
        doc="Units of the target catalog coordinates",
        default="degree",
    )


class AstropyMatchTask(pipeBase.Task):
    """A task for running the astropy matcher `match_to_catalog_sky` on
    between target and reference catalogs."""

    ConfigClass = AstropyMatchConfig

    def run(self, refCatalog, targetCatalog):
        """Run matcher

        Parameters
        ----------
        refCatalog: `pandas.core.frame.DataFrame`
            The reference catalog with coordinates in degrees
        targetCatalog: `pandas.core.frame.DataFrame`
            The target catalog with coordinates in degrees

        Returns
        -------
        `pipeBase.Struct` containing:
            refMatchIndices: `numpy.ndarray`
                Array of indices of matched reference catalog objects
            targetMatchIndices: `numpy.ndarray`
                Array of indices of matched target catalog objects
            separations: `astropy.coordinates.angles.Angle`
                Array of angle separations between matched objects
        """
        refCat_ap = SkyCoord(
            ra=refCatalog["coord_ra"].values * units.Unit(self.config.refCatUnits),
            dec=refCatalog["coord_dec"].values * units.Unit(self.config.refCatUnits),
        )

        sourceCat_ap = SkyCoord(
            ra=targetCatalog["coord_ra"].values * units.Unit(self.config.targetCatUnits),
            dec=targetCatalog["coord_dec"].values * units.Unit(self.config.targetCatUnits),
        )

        id, d2d, d3d = refCat_ap.match_to_catalog_sky(sourceCat_ap)

        goodMatches = d2d.arcsecond < self.config.maxDistance

        refMatchIndices = np.flatnonzero(goodMatches)
        targetMatchIndices = id[goodMatches]

        separations = d2d[goodMatches].arcsec

        return pipeBase.Struct(
            refMatchIndices=refMatchIndices, targetMatchIndices=targetMatchIndices, separations=separations
        )


class CatalogMatchConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"targetCatalog": "objectTable_tract", "refCatalog": "gaia_dr2_20200414"},
):
    catalog = pipeBase.connectionTypes.Input(
        doc="The tract-wide catalog to make plots from.",
        storageClass="DataFrame",
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
        storageClass="DataFrame",
        dimensions=("tract", "skymap"),
    )


class CatalogMatchConfig(pipeBase.PipelineTaskConfig, pipelineConnections=CatalogMatchConnections):
    matcher = pexConfig.ConfigurableField[pipeBase.Task](
        target=AstropyMatchTask, doc="Task for matching refCat and SourceCatalog"
    )

    epoch = pexConfig.Field[float](doc="Epoch to which reference objects are shifted", default=2015.0)

    bands = pexConfig.ListField[str](
        doc="All bands to persist to downstream tasks",
        default=["u", "g", "r", "i", "z", "y"],
    )

    selectorBand = pexConfig.Field[str](
        doc="Band to use when selecting objects, primarily for extendedness", default="i"
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
        default=["i_cModelFlux", "x", "y"],
    )

    requireProperMotion = pexConfig.Field[bool](
        doc="Only use reference catalog objects with proper motion information",
        default=False,
    )

    anyFilterMapsToThis = pexConfig.Field[str](
        doc="Any filter for the reference catalog maps to this",
        default="phot_g_mean",
    )


class CatalogMatchTask(pipeBase.PipelineTask):
    """Match a tract-level catalog to a reference catalog"""

    ConfigClass = CatalogMatchConfig
    _DefaultName = "analysisToolsCatalogMatch"

    def __init__(self, butler=None, initInputs=None, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("matcher")

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class

        inputs = butlerQC.get(inputRefs)

        columns = ["coord_ra", "coord_dec", "patch"] + self.config.extraColumns.list()
        for selectorAction in [
            self.config.selectorActions,
            self.config.sourceSelectorActions,
            self.config.extraColumnSelectors,
        ]:
            for selector in selectorAction:
                for band in self.config.bands:
                    selectorSchema = selector.getFormattedInputSchema(band=band)
                    columns += [s[0] for s in selectorSchema]

        dataFrame = inputs["catalog"].get(parameters={"columns": columns})
        inputs["catalog"] = dataFrame

        tract = butlerQC.quantum.dataId["tract"]

        self.refObjLoader = ReferenceObjectLoader(
            dataIds=[ref.datasetRef.dataId for ref in inputRefs.refCat],
            refCats=inputs.pop("refCat"),
            name=self.config.connections.refCat,
            log=self.log,
        )
        self.refObjLoader.config.requireProperMotion = self.config.requireProperMotion
        self.refObjLoader.config.anyFilterMapsToThis = self.config.anyFilterMapsToThis

        self.setRefCat(inputs.pop("skymap"), tract)

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    def run(self, catalog):
        """Prep the catalog and run the matcher.

        Parameters
        ----------
        catalog : `pandas.core.frame.DataFrame`

        Returns
        -------
        `pipeBase.Struct` containing:
            matchedCat : `pandas.core.frame.DataFrame`
                Catalog containing the matched objects with all columns from
                the original input catalogs, with the suffix "_ref" or
                "_target" for duplicated column names, plus a column with the
                angular separation in arcseconds between matches.
        """
        # Apply the selectors to the catalog
        mask = np.ones(len(catalog), dtype=bool)
        for selector in self.config.selectorActions:
            mask &= selector(catalog, bands=self.config.bands)

        for selector in self.config.sourceSelectorActions:
            mask &= selector(catalog, band=self.config.selectorBand).astype(bool)

        targetCatalog = catalog[mask]
        targetCatalog = targetCatalog.reset_index()

        if (len(targetCatalog) == 0) or (len(self.refCat) == 0):
            matches = pipeBase.Struct(
                refMatchIndices=np.array([]), targetMatchIndices=np.array([]), separations=np.array([])
            )
        else:
            # Run the matcher
            matches = self.matcher.run(self.refCat, targetCatalog)

        # Join the catalogs for the matched catalogs
        refMatches = self.refCat.iloc[matches.refMatchIndices].reset_index()
        sourceMatches = targetCatalog.iloc[matches.targetMatchIndices].reset_index()
        matchedCat = sourceMatches.join(refMatches, lsuffix="_target", rsuffix="_ref")

        separations = pd.Series(matches.separations).rename("separation")
        matchedCat = matchedCat.join(separations)

        return pipeBase.Struct(matchedCatalog=matchedCat)

    def setRefCat(self, skymap, tract):
        """Make a reference catalog with coordinates in degrees

        Parameters
        ----------
        skymap : `lsst.skymap`
            The skymap used to define the patch boundaries.
        tract : int
            The tract corresponding to the catalog data.
        """
        # Load the reference objects in a skyCircle around the tract
        tractInfo = skymap.generateTract(tract)
        boundingCircle = tractInfo.getOuterSkyPolygon().getBoundingCircle()
        center = lsst.geom.SpherePoint(boundingCircle.getCenter())
        radius = boundingCircle.getOpeningAngle()

        epoch = Time(self.config.epoch, format="decimalyear")

        skyCircle = self.refObjLoader.loadSkyCircle(center, radius, "i", epoch=epoch)
        refCat = skyCircle.refCat

        # Convert the coordinates to RA/Dec and convert the catalog to a
        # dataframe
        refCat["coord_ra"] = (refCat["coord_ra"] * units.radian).to(units.degree).to_value()
        refCat["coord_dec"] = (refCat["coord_dec"] * units.radian).to(units.degree).to_value()
        self.refCat = refCat.asAstropy().to_pandas()


class CatalogMatchVisitConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("visit",),
    defaultTemplates={"targetCatalog": "sourceTable_visit", "refCatalog": "gaia_dr2_20200414"},
):
    catalog = pipeBase.connectionTypes.Input(
        doc="The visit-wide catalog to make plots from.",
        storageClass="DataFrame",
        name="sourceTable_visit",
        dimensions=("visit",),
        deferLoad=True,
    )

    refCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The astrometry reference catalog to match to loaded input catalog sources.",
        name="{refCatalog}",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    visitSummaryTable = pipeBase.connectionTypes.Input(
        doc="A summary table of the ccds in the visit",
        storageClass="ExposureCatalog",
        name="finalVisitSummary",
        dimensions=("visit",),
    )

    matchedCatalog = pipeBase.connectionTypes.Output(
        doc="Catalog with matched target and reference objects with separations",
        name="{targetCatalog}_{refCatalog}_match",
        storageClass="DataFrame",
        dimensions=("visit",),
    )


class CatalogMatchVisitConfig(CatalogMatchConfig, pipelineConnections=CatalogMatchVisitConnections):
    selectorActions = ConfigurableActionStructField(
        doc="Which selectors to use to narrow down the data for QA plotting.",
        default={"flagSelector": VisitPlotFlagSelector()},
    )

    extraColumns = pexConfig.ListField[str](
        doc="Other catalog columns to persist to downstream tasks",
        default=["psfFlux", "psfFluxErr"],
    )

    def setDefaults(self):
        # sourceSelectorActions.sourceSelector is StarSelector
        self.sourceSelectorActions.sourceSelector.vectorKey = "extendedness"
        # extraColumnSelectors.selector1 is SnSelector
        self.extraColumnSelectors.selector1.fluxType = "psfFlux"
        # extraColumnSelectors.selector2 is GalaxySelector
        self.extraColumnSelectors.selector2.vectorKey = "extendedness"


class CatalogMatchVisitTask(CatalogMatchTask):
    """Match a visit-level catalog to a reference catalog"""

    ConfigClass = CatalogMatchVisitConfig
    _DefaultName = "analysisToolsCatalogMatchVisit"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class

        inputs = butlerQC.get(inputRefs)

        columns = ["coord_ra", "coord_dec", "detector"] + self.config.extraColumns.list()
        for selectorAction in [
            self.config.selectorActions,
            self.config.sourceSelectorActions,
            self.config.extraColumnSelectors,
        ]:
            for selector in selectorAction:
                selectorSchema = selector.getFormattedInputSchema()
                columns += [s[0] for s in selectorSchema]

        dataFrame = inputs["catalog"].get(parameters={"columns": columns})
        inputs["catalog"] = dataFrame

        self.refObjLoader = ReferenceObjectLoader(
            dataIds=[ref.datasetRef.dataId for ref in inputRefs.refCat],
            refCats=inputs.pop("refCat"),
            name=self.config.connections.refCat,
            log=self.log,
        )
        self.refObjLoader.config.requireProperMotion = self.config.requireProperMotion
        self.refObjLoader.config.anyFilterMapsToThis = self.config.anyFilterMapsToThis

        self.setRefCat(inputs.pop("visitSummaryTable"))

        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    def setRefCat(self, visitSummaryTable):
        """Make a reference catalog with coordinates in degrees

        Parameters
        ----------
        visitSummaryTable : `lsst.afw.table.ExposureCatalog`
            The table of visit information
        """
        # Get convex hull around the detectors, then get its center and radius
        corners = []
        for visSum in visitSummaryTable:
            for ra, dec in zip(visSum["raCorners"], visSum["decCorners"]):
                corners.append(lsst.geom.SpherePoint(ra, dec, units=lsst.geom.degrees).getVector())
        visitBoundingCircle = lsst.sphgeom.ConvexPolygon.convexHull(corners).getBoundingCircle()
        center = lsst.geom.SpherePoint(visitBoundingCircle.getCenter())
        radius = visitBoundingCircle.getOpeningAngle()

        # Get the observation date of the visit
        obsDate = visSum.getVisitInfo().getDate()
        epoch = Time(obsDate.toPython())

        # Load the reference catalog in the skyCircle of the detectors, then
        # convert the coordinates to degrees and convert the catalog to a
        # dataframe
        skyCircle = self.refObjLoader.loadSkyCircle(center, radius, "i", epoch=epoch)
        refCat = skyCircle.refCat

        refCat["coord_ra"] = (refCat["coord_ra"] * units.radian).to(units.degree).to_value()
        refCat["coord_dec"] = (refCat["coord_dec"] * units.radian).to(units.degree).to_value()
        self.refCat = refCat.asAstropy().to_pandas()
