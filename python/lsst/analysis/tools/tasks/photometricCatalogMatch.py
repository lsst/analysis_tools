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

__all__ = ("PhotometricCatalogMatchConfig", "PhotometricCatalogMatchTask")


import numpy as np
from astropy.time import Time
from astropy.table import Table, hstack
from smatch import Matcher

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask
from lsst.skymap import BaseSkyMap
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
import lsst.geom
from ..actions.vector import (
    CoaddPlotFlagSelector,
    VisitPlotFlagSelector,
    GalaxySelector,
    SnSelector,
    StarSelector,
)

from ..tasks.catalogMatch import CatalogMatchConfig, CatalogMatchTask, CatalogMatchConnections

class PhotometricCatalogMatchConnections(CatalogMatchConnections):
    pass

class PhotometricCatalogMatchConfig(
    CatalogMatchConfig, pipelineConnections=PhotometricCatalogMatchConnections
):

    filterNames = pexConfig.ListField[str](
        doc="Physical filter names to persist downstream.",
        default=["HSC-G", "HSC-R", "HSC-I", "HSC-Z", "HSC-Y"],
    )

    # in the validation, we need to confirm that all of the filter names are in the filter map.

    # do a setDefaults and set however you like.
    def setDefaults(self):
        super().setDefaults()
        self.referenceCatalogLoader.doReferenceSelection = False
        self.referenceCatalogLoader.doApplyColorTerms = True


class PhotometricCatalogMatchTask(CatalogMatchTask):
    """A wrapper task to provide the information that
    is specific to the photometric reference catalog.
    """
    ConfigClass = PhotometricCatalogMatchConfig
    _DefaultName = "analysisToolsPhotometricCatalogMatch"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Run the matching to the photometric reference
        catalog.

        Parameters
        ----------
        `butlerQC` : lsst.pipe.base.butlerQuantumContext.ButlerQuantumContext
        `inputRefs` : lsst.pipe.base.connections.InputQuantizedConnection
        `outputRefs` : lsst.pipe.base.connections.OutputQuantizedConnection

        """

        inputs = butlerQC.get(inputRefs)
        bands = []
        for filterName in self.config.filterNames:
            bands.append(self.config.referenceCatalogLoader.refObjLoader.filterMap[filterName])

        columns = self.prepColumns(bands)
        table = inputs["catalog"].get(parameters={"columns": columns})
        inputs["catalog"] = table

        tract = butlerQC.quantum.dataId["tract"]

        loaderTask = LoadReferenceCatalogTask(
            config=self.config.referenceCatalogLoader,
            dataIds=[ref.dataId for ref in inputRefs.refCat],
            name=inputs["refCat"][0].ref.datasetType.name,
            refCats=inputs["refCat"],
        )

        skymap = inputs.pop("skymap")
        loadedRefCat = self._loadRefCat(loaderTask, skymap[tract])
        outputs = self.run(catalog=inputs["catalog"], loadedRefCat=loadedRefCat, bands=bands)

        butlerQC.put(outputs, outputRefs)


class PhotometricCatalogMatchVisitConnections(
    CatalogMatchConnections,
    dimensions=("visit",),
    defaultTemplates={"targetCatalog": "sourceTable_visit", "refCatalog": "ps1_pv3_3pi_20170110"},
    ):
    
    catalog = pipeBase.connectionTypes.Input(
        doc="The visit-wide catalog to make plots from.",
        storageClass="ArrowAstropy",
        name="{targetCatalog}",
        dimensions=("visit",),
        deferLoad=True,
    )

    refCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The photometric reference catalog to match to.",
        name="{refCatalog}",
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



class PhotometricCatalogMatchVisitConfig(
    PhotometricCatalogMatchConfig, pipelineConnections=PhotometricCatalogMatchVisitConnections
):

    def setDefaults(self):
        self.filterNames = []
        self.extraPerBandColumns = []
        self.patchColumn = ""
        self.selectorBands = []
        self.selectorActions.flagSelector = VisitPlotFlagSelector
        self.sourceSelectorActions.sourceSelector.vectorKey = "extendedness"
        self.extraColumnSelectors.selector1.fluxType = "psfFlux"
        self.extraColumnSelectors.selector2.vectorKey = "extendedness"

class PhotometricCatalogMatchVisitTask(PhotometricCatalogMatchTask):
    """A wrapper task to provide the information that
    is specific to the photometric reference catalog.
    """
    ConfigClass = PhotometricCatalogMatchVisitConfig
    _DefaultName = "analysisToolsPhotometricCatalogMatchVisit"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Run the matching to the photometric reference
        catalog.

        Parameters
        ----------
        `butlerQC` : lsst.pipe.base.butlerQuantumContext.ButlerQuantumContext
        `inputRefs` : lsst.pipe.base.connections.InputQuantizedConnection
        `outputRefs` : lsst.pipe.base.connections.OutputQuantizedConnection

        """

        inputs = butlerQC.get(inputRefs)
        physicalFilter = inputs["catalog"].dataId["physical_filter"]
        bands = [self.config.referenceCatalogLoader.refObjLoader.filterMap[physicalFilter]]
        # No bands needed for visit tables
        # but we do need them later for the matching
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
        inputs["catalog"] = table

        loaderTask = LoadReferenceCatalogTask(
            config=self.config.referenceCatalogLoader,
            dataIds=[ref.dataId for ref in inputRefs.refCat],
            name=inputs["refCat"][0].ref.datasetType.name,
            refCats=inputs["refCat"],
        )

        visitSummaryTable = inputs.pop("visitSummaryTable")
        loadedRefCat = self._loadRefCat(loaderTask, visitSummaryTable, physicalFilter)
        outputs = self.run(catalog=inputs["catalog"], loadedRefCat=loadedRefCat, bands=bands)

        # The matcher adds the band to the front of the columns
        # but the visit plots aren't expecting it
        matchedCat = outputs.matchedCatalog
        cols = list(outputs.matchedCatalog.columns)
        for col in cols:
            if col[:2] == bands[0] + "_":
                outputs.matchedCatalog.rename_column(col, col[2:])

        butlerQC.put(outputs, outputRefs)

    def _loadRefCat(self, loaderTask, visitSummaryTable, physicalFilter):
        """Make a reference catalog with coordinates in degrees

        Parameters
        ----------
        visitSummaryTable : `lsst.afw.table.ExposureCatalog`
            The table of visit information
        """
        # Get convex hull around the detectors, then get its center and radius
        corners = []
        for visSum in visitSummaryTable:
            for (ra, dec) in zip(visSum["raCorners"], visSum["decCorners"]):
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

        loadedRefCat = loaderTask.getSkyCircleCatalog(center, radius, [physicalFilter], epoch=epoch)

        return Table(loadedRefCat)
