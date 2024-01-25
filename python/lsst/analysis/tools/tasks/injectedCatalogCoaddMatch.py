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

__all__ = ("InjectedCatalogCoaddMatchConfig", "InjectedCatalogCoaddMatchTask")

import lsst.pipe.base as pipeBase
import numpy as np
from lsst.pex.config import ListField
from lsst.pex.config.configurableActions import ConfigurableActionStructField
from lsst.pipe.base.connectionTypes import Input, Output

from ..actions.vector.selectors import CoaddPlotFlagSelector, InjectedCoaddPlotFlagSelector, StarSelector
from ..interfaces import VectorAction
from ..tasks.catalogMatch import CatalogMatchConnections, CatalogMatchConfig, CatalogMatchTask


# TODO: Inherit from CatalogMatchConnections
class InjectedCatalogCoaddMatchConnections(
    CatalogMatchConnections,
    # dimensions=("tract", "skymap"),
    # defaultTemplates={
    #     "injectedTargetCatalog": "injected_objectTable_tract",
    #     "injectedReferenceCatalog": "injected_deepCoadd_catalog_tract",
    # },
):
    # targetCatalog = Input(
    #     doc="The tract-wide catalog to make plots from.",
    #     storageClass="ArrowAstropy",
    #     name="{injectedTargetCatalog}",
    #     dimensions=("tract", "skymap"),
    #     deferLoad=True,
    # )
    refCat = Input(
        doc="Tract level catalog of injected sources.",
        name="{refCatalog}",
        storageClass="ArrowAstropy",
        dimensions=("skymap", "tract"),
    )

    # matchedCatalog = Output(
    #     doc="Catalog with matched target and reference objects with separations",
    #     name="{injectedTargetCatalog}_{injectedReferenceCatalog}_match",
    #     storageClass="ArrowAstropy",
    #     dimensions=("tract", "skymap"),
    # )
    # def __init__(self, *, config=None):
    #     super().__init__(config=config)
    #     self.inputs.remove("refCat")


class InjectedCatalogCoaddMatchConfig(
    CatalogMatchConfig, pipelineConnections=InjectedCatalogCoaddMatchConnections
):
    bands = ListField[str](
        doc="List of photometric bands used in the reference and target catalogs.",
        default=["g", "r", "i", "z", "y"],
    )
    sourceSelectorActions = ConfigurableActionStructField[VectorAction](
        doc="Which selectors to use to narrow down the data for QA plotting.",
        default={
            "sourceSelector1": CoaddPlotFlagSelector(),
            "sourceSelector2": StarSelector(),
        },
    )
    injectedSourceSelectorActions = ConfigurableActionStructField[VectorAction](
        doc="What types of injected sources to use.",
        default={
            "injectedSourceSelector1": InjectedCoaddPlotFlagSelector(),
        },
    )

    def setDefaults(self):
        super().setDefaults()
        self.connections.targetCatalog = "injected_objectTable_tract"
        self.connections.refCatalog = "injected_deepCoadd_catalog_tract"
        self.returnNonMatches = True


class InjectedCatalogCoaddMatchTask(CatalogMatchTask):
    """A wrapper task to provide the information that
    is specific to the photometric reference catalog.
    """

    ConfigClass = InjectedCatalogCoaddMatchConfig
    _DefaultName = "analysisToolsPhotometricCatalogMatch"

    # TODO: make definition for flag selection on the refCat

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

        columns = self.prepColumns(self.config.bands)
        targetCat = inputs["catalog"].get(parameters={"columns": columns})
        injectedCat = inputs["refCat"]
        # Apply the selectors to the injected catalog
        mask = np.ones(len(injectedCat), dtype=bool)
        for selector in self.config.injectedSourceSelectorActions:
            mask &= selector(injectedCat).astype(bool)

        injectedCat = injectedCat[mask]
        outputs = self.run(targetCatalog=targetCat, refCatalog=injectedCat, bands=self.config.bands)

        butlerQC.put(outputs, outputRefs)
