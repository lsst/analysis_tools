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

__all__ = (
    "AstrometricCatalogMatchConfig",
    "AstrometricCatalogMatchTask",
    "AstrometricCatalogMatchVisitConfig",
    "AstrometricCatalogMatchVisitTask",
)

import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import numpy as np
from astropy.table import Table
from astropy.time import Time
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask

from ..actions.vector import VisitPlotFlagSelector
from ..tasks.catalogMatch import CatalogMatchConfig, CatalogMatchConnections, CatalogMatchTask


class AstrometricCatalogMatchConfig(CatalogMatchConfig, pipelineConnections=CatalogMatchConnections):
    bands = pexConfig.ListField[str](
        doc="The bands to persist downstream",
        default=["u", "g", "r", "i", "z", "y"],
    )

    def setDefaults(self):
        super().setDefaults()
        self.referenceCatalogLoader.doApplyColorTerms = False
        self.referenceCatalogLoader.refObjLoader.requireProperMotion = True
        self.referenceCatalogLoader.refObjLoader.anyFilterMapsToThis = "phot_g_mean"


class AstrometricCatalogMatchTask(CatalogMatchTask):
    """Match a tract-level catalog to a reference catalog"""

    ConfigClass = AstrometricCatalogMatchConfig
    _DefaultName = "analysisToolsAstrometricCatalogMatch"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class

        inputs = butlerQC.get(inputRefs)
        columns = self.prepColumns(self.config.bands)
        table = inputs["catalog"].get(parameters={"columns": columns})

        tract = butlerQC.quantum.dataId["tract"]

        loaderTask = LoadReferenceCatalogTask(
            config=self.config.referenceCatalogLoader,
            dataIds=[ref.dataId for ref in inputRefs.refCat],
            name=inputs["refCat"][0].ref.datasetType.name,
            refCats=inputs["refCat"],
        )

        skymap = inputs.pop("skymap")
        loadedRefCat = self._loadRefCat(loaderTask, skymap[tract])
        outputs = self.run(catalog=table, loadedRefCat=loadedRefCat, bands=self.config.bands)

        butlerQC.put(outputs, outputRefs)


class AstrometricCatalogMatchVisitConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("visit",),
    defaultTemplates={"targetCatalog": "sourceTable_visit", "refCatalog": "gaia_dr2_20200414"},
):
    catalog = pipeBase.connectionTypes.Input(
        doc="The visit-wide catalog to make plots from.",
        storageClass="ArrowAstropy",
        name="sourceTable_visit",
        dimensions=("visit",),
        deferLoad=True,
    )

    refCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The astrometry reference catalog to match to loaded input catalog sources.",
        name="gaia_dr2_20200414",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    visitSummaryTable = pipeBase.connectionTypes.Input(
        doc="A summary table of the ccds in the visit",
        storageClass="ExposureCatalog",
        name="visitSummary",
        dimensions=("visit",),
    )

    matchedCatalog = pipeBase.connectionTypes.Output(
        doc="Catalog with matched target and reference objects with separations",
        name="{targetCatalog}_{refCatalog}_match",
        storageClass="ArrowAstropy",
        dimensions=("visit",),
    )


class AstrometricCatalogMatchVisitConfig(
    AstrometricCatalogMatchConfig, pipelineConnections=AstrometricCatalogMatchVisitConnections
):
    selectorActions = ConfigurableActionStructField(
        doc="Which selectors to use to narrow down the data for QA plotting.",
        default={"flagSelector": VisitPlotFlagSelector},
    )

    extraColumns = pexConfig.ListField[str](
        doc="Other catalog columns to persist to downstream tasks",
        default=["psfFlux", "psfFluxErr"],
    )

    bands = pexConfig.ListField[str](
        doc="The bands to persist downstream",
        default=[],
    )

    def setDefaults(self):
        # sourceSelectorActions.sourceSelector is StarSelector
        self.sourceSelectorActions.sourceSelector.vectorKey = "extendedness"
        # extraColumnSelectors.selector1 is SnSelector
        self.extraColumnSelectors.selector1.fluxType = "psfFlux"
        # extraColumnSelectors.selector2 is GalaxySelector
        self.extraColumnSelectors.selector2.vectorKey = "extendedness"
        self.referenceCatalogLoader.doApplyColorTerms = False
        self.referenceCatalogLoader.refObjLoader.requireProperMotion = False
        self.referenceCatalogLoader.refObjLoader.anyFilterMapsToThis = "phot_g_mean"


class AstrometricCatalogMatchVisitTask(AstrometricCatalogMatchTask):
    """Match a visit-level catalog to a reference catalog"""

    ConfigClass = AstrometricCatalogMatchVisitConfig
    _DefaultName = "analysisToolsAstrometricCatalogMatchVisit"

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

        table = inputs["catalog"].get(parameters={"columns": columns})

        loaderTask = LoadReferenceCatalogTask(
            config=self.config.referenceCatalogLoader,
            dataIds=[ref.dataId for ref in inputRefs.refCat],
            name=inputs["refCat"][0].ref.datasetType.name,
            refCats=inputs["refCat"],
        )

        visitSummaryTable = inputs.pop("visitSummaryTable")
        loadedRefCat = self._loadRefCat(loaderTask, visitSummaryTable)
        outputs = self.run(catalog=table, loadedRefCat=loadedRefCat, bands=self.config.bands)

        butlerQC.put(outputs, outputRefs)

    def _loadRefCat(self, loaderTask, visitSummaryTable):
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
                # If the coordinates are nan then don't keep going
                # because it crashes later
                if not np.isfinite(ra) or not np.isfinite(dec):
                    raise pipeBase.NoWorkFound("Visit summary corners not finite")
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

        filterName = self.config.referenceCatalogLoader.refObjLoader.anyFilterMapsToThis
        loadedRefCat = loaderTask.getSkyCircleCatalog(center, radius, filterName, epoch=epoch)

        return Table(loadedRefCat)
