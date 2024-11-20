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

__all__ = ("AssociatedSourcesCatalogMatchConfig", "AssociatedSourcesCatalogMatchTask")


import lsst.geom
import lsst.pipe.base as pipeBase
import numpy as np
from astropy.table import Table
from astropy.time import Time
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask

from ..actions.vector import StarSelector, VisitPlotFlagSelector
from ..tasks.catalogMatch import CatalogMatchConfig, CatalogMatchConnections, CatalogMatchTask


class AssociatedSourcesCatalogMatchConnections(
    CatalogMatchConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"targetCatalog": "isolated_star_presources", "refCatalog": "the_monster_20240904"},
):
    matchedCatalog = pipeBase.connectionTypes.Output(
        doc="Catalog with matched target and reference objects with separations",
        name="{targetCatalog}_{refCatalog}_match_assoc",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )
class AssociatedSourcesCatalogMatchConfig(
    CatalogMatchConfig,
    pipelineConnections=AssociatedSourcesCatalogMatchConnections,
):
    def setDefaults(self):
        super().setDefaults()
        self.matchesRefCat = True
        self.referenceCatalogLoader.doReferenceSelection = False
        self.referenceCatalogLoader.doApplyColorTerms = False


class AssociatedSourcesCatalogMatchTask(CatalogMatchTask):
    """A wrapper task to provide the information that
    is specific to the associated sources reference catalog.
    """

    ConfigClass = AssociatedSourcesCatalogMatchConfig
    _DefaultName = "analysisToolsAssociatedSourcesCatalogMatch"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Run the matching to the photometric reference catalog.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
        inputRefs : `lsst.pipe.base.InputQuantizedConnection`
        outputRefs : `lsst.pipe.base.OutputQuantizedConnection`
        """

        inputs = butlerQC.get(inputRefs)
        bands = []
        for filterName in self.config.filterNames:
            bands.append(self.config.referenceCatalogLoader.refObjLoader.filterMap[filterName])

        # For some reason the imsim filterMaps don't work the same way
        # as the HSC ones do, this is a bit hacky but fixes this.
        # This code makes the assumption that filterMap is a dict
        # mapping observed filter names: band, but filterMap is
        # currently defined as a mapping of:
        # observed filter names -> reference catalog filter names.
        # Therefore, when the reference catalog filter name is
        # not the band name we need this work around.
        # TODO: workaround for DM-46728
        if bands[0].startswith("lsst") or "sim" in bands[0] or "smeared" in bands[0]:
            bands = self.config.filterNames
        elif bands[0].startswith("monster"):
            # for the_monster_20240904 the reference catalog filter name is
            # "monster_{system}_{band}" the last character is the band
            bands = [band[-1] for band in bands]

        columns = self.prepColumns(bands)
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
        # subset cat to ra/dec of unique sources
        uniqueObj, uniqueInd = np.unique(table['obj_index'], return_index=True)
        uniqueTable = table[uniqueInd].copy()
        outputs = self.run(targetCatalog=uniqueTable, refCatalog=loadedRefCat, bands=bands)
        import pdb; pdb.set_trace()
        #make a new empty table with columns of outputs['matchedCatalog']
        newColsTable = Table({col: [] for col in outputs['matchedCatalog'].colnames})
        for obj in uniqueObj:
            objInds = np.where(table['object_index'] == obj)[0]

            # objTable = table[objInd].copy()
            # objOutputs = self.run(targetCatalog=objTable, refCatalog=loadedRefCat, bands=bands)
            # outputs['matchedCatalog'] = outputs['matchedCatalog'].copy().append(objOutputs['matchedCatalog'])
        outputs['matchedCatalog']
        # expand out subset to full catalog and add to outputs

        butlerQC.put(outputs, outputRefs)