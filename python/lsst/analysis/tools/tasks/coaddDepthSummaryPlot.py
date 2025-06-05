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
from astropy.table import Table

import numpy as np

__all__ = (
    "CoaddDepthSummaryPlotConfig",
    "CoaddDepthSummaryPlotTask",
)
from lsst.pex.config import ListField
from lsst.pipe.base import connectionTypes as cT
from lsst.pipe.base import Struct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CoaddDepthSummaryPlotConnections(
    AnalysisBaseConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"inputName": "coadd_depth_table",
                      "outputName": "coadd_depth_per_pixel",},
):
    statTable = cT.Input(
        doc="Table with resulting n_image based depth statistics.",
        name="{inputName}",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )
    
    data = cT.Input(
        doc="Coadd n_image to load from the butler (pixel values are the number of input images).",
        name="template_coadd_n_image",
        storageClass="ImageU",
        multiple=True,
        dimensions=("tract", "patch", "band", "skymap"),
        deferLoad=True,
    )
    
    depth_table = cT.Output(
        doc="Table with coadd depth per pixel.",
        name="{outputName}",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )
    
    calculated_quantiles = cT.Output(
        doc="Table with calculated percentile values.",
        name="calculated_quantiles",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )


class CoaddDepthSummaryPlotConfig(AnalysisBaseConfig,
                                  pipelineConnections=CoaddDepthSummaryPlotConnections):
    quantile_list = ListField(
        default=[5, 10, 25, 50, 75, 90, 95],
        dtype=int,
        doc="The percentiles at which to compute n_image values, in ascending order.",
    )


class CoaddDepthSummaryPlotTask(AnalysisPipelineTask):
    ConfigClass = CoaddDepthSummaryPlotConfig
    _DefaultName = "coaddDepthSummaryPlot"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputs):
        bands = []
        patches = []
        
        q_bands = []
        q_patches = []
        
        depths = []
        pixels = []
        quantiles = []

        for n_image_handle in inputs["data"]:
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
            
            quantile = list(np.percentile(n_image.array, q=self.config.quantile_list))

            quantiles.extend(quantile)

            for i in range(len(quantile)):
                q_bands.append(band)
                q_patches.append(patch)

        names = ["patch", "band", "depth", "pixels"]
        data = [patches, bands, depths, pixels]

        q_names = ["q_patch", "q_band", "quantiles"]
        q_data = [q_patches, q_bands, quantiles]
        
        depth_table = Table(data=data, names=names)
        calculated_quantiles = Table(data=q_data, names=q_names)
        
        return Struct(depth_table=depth_table, calculated_quantiles=calculated_quantiles)
