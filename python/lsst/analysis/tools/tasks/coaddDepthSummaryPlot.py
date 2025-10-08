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
    "CoaddDepthSummaryPlotConfig",
    "CoaddDepthSummaryPlotTask",
)

import numpy as np
from lsst.pipe.base import Struct
from lsst.pipe.base import connectionTypes as cT
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask, KeyedData


class CoaddDepthSummaryPlotConnections(
    AnalysisBaseConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={
        "coaddName": "",
    },
):
    skymap = cT.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    n_image_data = cT.Input(
        doc="Coadd n_image to load from the butler (pixel values are the number of input images).",
        name="{coaddName}_coadd_n_image",
        storageClass="ImageU",
        multiple=True,
        dimensions=("tract", "patch", "band", "skymap"),
        deferLoad=True,
    )


class CoaddDepthSummaryPlotConfig(AnalysisBaseConfig, pipelineConnections=CoaddDepthSummaryPlotConnections):
    pass


class CoaddDepthSummaryPlotTask(AnalysisPipelineTask):
    ConfigClass = CoaddDepthSummaryPlotConfig
    _DefaultName = "coaddDepthSummaryPlot"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        skymap = inputs["skymap"]
        dataId = butlerQC.quantum.dataId
        tractInfo = skymap[dataId["tract"]]
        outputs = self.run(data={"n_image_data": inputs["n_image_data"]}, tractInfo=tractInfo)
        butlerQC.put(outputs, outputRefs)

    def run(self, *, data: KeyedData | None = None, **kwargs) -> Struct:
        """Use n_images to make a plot illustrating coadd depth."""
        bands = []
        patches = []

        depths = []
        pixels = []

        for n_image_handle in data["n_image_data"]:
            n_image = n_image_handle.get()
            data_id = n_image_handle.dataId
            band = str(data_id.band.name)
            patch = int(data_id.patch.id)

            depth, pixel = np.unique(n_image.array, return_counts=True)

            depths.extend(depth)
            pixels.extend(pixel)

            for i in range(len(depth)):
                bands.append(band)
                patches.append(patch)

        pixel_data = {"patch": patches, "band": bands, "depth": depths, "pixels": pixels}

        outputs = super().run(data=pixel_data, **kwargs)  # this creates a struct for the output

        return outputs
