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
from smatch import Matcher

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask


class PhotometricCatalogMatchConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"targetCatalog": "objectTable_tract",
                      "refCatalog": "ps1_whatever"},
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
        name="ps1_whatever",
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


class PhotometricCatalogMatchConfig(
    pipeBase.PipelineTaskConfig,
    pipelineConnections=PhotometricCatalogMatchConnections
):

    referenceCatalogLoader = pexConfig.ConfigurableField(
        target=LoadReferenceCatalogTask, doc="Reference catalog loader",
    )

    epoch = pexConfig.Field[float](
        doc="Epoch to which reference objects are shifted.",
        default=2015.0,
    )

    # bands = pexConfig.ListField[str](
    #     doc="All bands to persist to downstream tasks.",
    #     default=["u", "g", "r", "i", "z", "y"],
    # )
    filterNames = pexConfig.ListField[str](
        doc="Physical filter names to persist downstream.",
        default=["HSC-G", "HSC-R", "HSC-I", "HSC-Z", "HSC-Z"],
    )

    selectorBand = pexConfig.Field[str](
        doc="Band to use when selecting objects, primarily for extendedness.",
        default="i",
    )

    selectorActions = ConfigurableActionStructField(
        doc="Which selectors to use to narrow down the data for QA plotting.",
        default={"flagSelector": CoaddPlotFlagSelector},
    )

    sourceSelectorActions = ConfigurableActionStructField(
        doc="What types of sources to use.",
        default={"sourceSelector": StarSelector},
    )

    extraColumnSelectors = ConfigurableActionStructField(
        doc="Other selectors that are not used in this task, but whose columns; may be needed downstream.",
        default={"selector1": SnSelector, "selector2": GalaxySelector},
    )

    extraColumns = pexConfig.ListField[str](
        doc="Other catalog columns to persist to downstream tasks.",
        default=["i_cModelFlux", "x", "y"],
    )

    raColumn = pexConfig.Field[str](doc="RA column.", default="coord_ra")
    decColumn = pexConfig.Field[str](doc="Dec column.", default="coord_dec")
    patchColumn = pexConfig.Field[str](doc="Patch column.", default="patch")

    # in the validation, we need to confirm that all of the filter names are in the filter map.

    # do a setDefaults and set however you like.


class PhotometricCatalogMatchTask(pipeBase.PipelineTask):
    ConfigClass = PhotometricCatalogMatchConfig
    _DefaultName = "analysisToolsPhotometricCatalogMatch"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Make any subtasks?  matching?

    def runQuantum(self, butlerQC, inputRefs, outputRefs):

        inputs = butlerQC.get(inputRefs)

        columns = [self.config.raColumn, self.config.decColumn, self.config.patchColumn] + self.config.extraColumns.list()
        for selectorAction in [
            self.config.selectorActions,
            self.config.sourceSelectorActions,
            self.config.extraColumnSelectors,
        ]:
            for selector in selectorAction:
                for band in self.config.bands:
                    selectorSchema = selector.getFormattedInputSchema(band=band)
                    columns += [s[0] for s in selectorSchema]

        table = inputs["catalog"].get(parameters={"columns": columns})
        inputs["catalog"] = table

        tract = butlerQC.quantum.dataId["tract"]

        loaderTask = LoadReferenceCatalogTask(
            config=self.config.referenceCatalogLoader,
            dataIds=[ref.datasetRef for ref in inputRefs.refCat],
            refCats=inputs.pop("refCat"),
            name=self.config.connections.refCat,
        )

        skymap = inputs.pop("skymap")
        loadedRefCat = self._loadRefCat(loaderTask, skymap[tract])

        outputs = self.run(catalog=inputs["catalog"], loadedRefCat=loadedRefCat)

        butlerQC.put(outputs, outputRefs)

    def _loadRefCat(self, loaderTask, tractInfo):
        """docstring"""
        boundingCircle = tractInfo.getOuterSkyPolygon().getBoundingCircle()
        center = lsst.geom.SpherePoint(boundingCircle.getCenter())
        radius = boundingCircle.getOpeningAngle()

        epoch = Time(self.config.epoch, format="decimalyear")

        # This is always going to return degrees.
        loadedRefCat = loaderTask.getSkyCircleCatalog(center, radius, self.config.filterNames, epoch=epoch)

        return loadedRefCat

    def run(self, *, catalog, loadedRefCat):
        """docstring"""
        # Apply the selectors to the catalog
        mask = np.ones(len(catalog), dtype=bool)
        for selector in self.config.selectorActions:
            mask &= selector(catalog, bands=self.config.bands)

        for selector in self.config.sourceSelectorActions:
            mask &= selector(catalog, band=self.config.selectorBand).astype(bool)

        targetCatalog = catalog[mask]

        if (len(targetCatalog) == 0) or (len(loadedRefCat) == 0):
            # matches = pipeBase.Strict(
            #     refMatchIndices=np.array([]), targetMatchIndices=np.array([]), separations=np.array([])
            # )
            refMatchIndices = np.array([], dtype=np.int64)
            targetMatchIndices = np.array([], dtype=np.int64)
            dist = np.array([], dtype=np.float64)
        else:
            # Run the matcher.

            # This all assumes that everything is in degrees.
            # Which I think is okay, but the current task allows different units.
            # Need to configure match radius, either in this task or a subtask.
            with Matcher(targetCatalog["ra"], targetCatalog["dec"]) as m:
                idx, targetMatchIndices, refMatchIndices, dist = m.query_radius(loadedRefCat["ra"], loadedRefCat["dec"], 1./3600.)
            # Convert degrees to arcseconds.
            dist *= 3600.

        # And here we crop the catalog to the matches, and add columns (hstack?) with the reference mags.
        # And the reference positions (with new names).
        # Rename everything.
        # loadedRefCat["refMag"] has shape (Nstar, Nfilter)
        # loadedRefCat["refMag"] has same shape.
        # Can separate this out with band names which come from filterMap.
        # self.config.referenceCatalogLoader.refObjLoader.filterMap  That's fun!

        # Make sure you add the dist column!

        # return a struct with the output catalog.

        # obs_subaru/config/analysisToolsPhotometricCatalogMatch.py
        # config.referenceCatalogLoader.applyColorTerms = True
        # config.referenceCatalogLoader.colorterms.load(os.path.join(OBS_CONFIG_DIR, "colorterms.py"))
        # config.referenceCatalogLoader.refObjLoader.load(os.path.join(OBS_CONFIG_DIR, "filterMap.py"))
