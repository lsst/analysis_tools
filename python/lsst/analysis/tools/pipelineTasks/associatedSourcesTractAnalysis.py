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

from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

# These need to be updated for this analysis context
# from ..analysisPlots.analysisPlots import ShapeSizeFractionalDiffScatter
# from ..analysisPlots.analysisPlots import Ap12_PSF_skyPlot
# from ..analysisMetrics.analysisMetrics import ShapeSizeFractionalMetric


class AssociatedSourcesTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={
        "inputName": "isolated_star_sources",
        # "associatedSourcesInputName": "isolated_star_sources"},
    },
):
    data = ct.Input(
        doc="Visit based source table to load from the butler",
        name="sourceTable_visit",
        storageClass="DataFrame",
        # deferLoad=True,
        dimensions=("visit", "band"),
        multiple=True,
    )

    associatedSources = ct.Input(
        doc="Table of associated sources",
        # name="{associatedSourcesInputName}",
        name="{inputName}",
        storageClass="DataFrame",
        # deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
    )

    skyMap = ct.Input(
        doc="Input definition of geometry/bbox and projection/wcs for warped exposures",
        name="skyMap",
        storageClass="SkyMap",
        dimensions=("skymap",),
    )


class AssociatedSourcesTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssociatedSourcesTractAnalysisConnections
):
    def setDefaults(self):
        super().setDefaults()
        # set plots to run
        # update for this analysis context
        # self.plots.shapeSizeFractionalDiffScatter = \
        #     ShapeSizeFractionalDiffScatter()
        # self.plots.Ap12_PSF_skyPlot = Ap12_PSF_skyPlot()

        # set metrics to run
        # update for this analysis context
        # self.metrics.shapeSizeFractionalMetric = ShapeSizeFractionalMetric()


class AssociatedSourcesTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AssociatedSourcesTractAnalysisConfig
    _DefaultName = "associatedSourcesTractAnalysisTask"

    def getBoxWcs(self, skymap, tract):
        tractInfo = skymap.generateTract(tract)
        wcs = tractInfo.getWcs()
        tractBox = tractInfo.getBBox()
        self.log.info("Running tract: %s", tract)
        return tractBox, wcs

    def prepareAssociatedSources(self, skymap, tract, sourceCatalogs, associatedSources):
        """
        This should be a standalone function rather than being associated with
        this class.
        """

        import lsst.geom as geom
        import numpy as np
        import pandas as pd

        # Keep only sources with associations
        dataJoined = pd.concat(sourceCatalogs).merge(associatedSources, on="sourceId", how="inner")
        dataJoined.set_index("sourceId", inplace=True)

        # Determine which sources are contained in tract
        ra = np.radians(dataJoined["coord_ra"].values)
        dec = np.radians(dataJoined["coord_dec"].values)
        box, wcs = self.getBoxWcs(skymap, tract)
        box = geom.Box2D(box)
        x, y = wcs.skyToPixelArray(ra, dec)
        boxSelection = box.contains(x, y)

        # Keep only the sources in groups that are fully contained within the
        # tract
        dataJoined["boxSelection"] = boxSelection
        dataFiltered = dataJoined.groupby("obj_index").filter(lambda x: all(x["boxSelection"]))
        dataFiltered.drop(columns="boxSelection", inplace=True)

        return dataFiltered

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dataFiltered = self.prepareAssociatedSources(
            inputs["skyMap"],
            inputRefs.associatedSources.dataId.byName()["tract"],
            inputs["data"],
            inputs["associatedSources"],
        )

        kwargs = {"data": dataFiltered}

        outputs = self.run(**kwargs)
        butlerQC.put(outputs, outputRefs)
