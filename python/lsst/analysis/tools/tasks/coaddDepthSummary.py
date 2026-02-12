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
    "CoaddDepthSummaryConfig",
    "CoaddDepthSummaryTask",
)


import numpy as np
from astropy.table import Table

from lsst.pex.config import ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base import connectionTypes as cT


class CoaddDepthSummaryConnections(
    PipelineTaskConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"coaddName": ""},  # set as either deep or template in the pipeline
):
    data = cT.Input(
        doc="Coadd n_image to load from the butler (pixel values are the number of input images).",
        name="{coaddName}_coadd_n_image",
        storageClass="ImageU",
        multiple=True,
        dimensions=("tract", "patch", "band", "skymap"),
        deferLoad=True,
    )

    statTable = cT.Output(
        doc="Table with resulting n_image based depth statistics.",
        name="{coaddName}_coadd_depth_table",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )


class CoaddDepthSummaryConfig(PipelineTaskConfig, pipelineConnections=CoaddDepthSummaryConnections):
    threshold_list = ListField(
        default=[1, 3, 5, 12],
        dtype=int,
        doc="The n_image pixel value thresholds, in ascending order.",
    )

    quantile_list = ListField(
        default=[5, 10, 25, 50, 75, 90, 95],
        dtype=int,
        doc="The percentiles at which to compute n_image values, in ascending order.",
    )


class CoaddDepthSummaryTask(PipelineTask):
    ConfigClass = CoaddDepthSummaryConfig
    _DefaultName = "coaddDepthSummary"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, inputs):
        t = Table()
        bands = []
        patches = []
        medians = []
        stdevs = []
        stats = []
        quantiles = []

        for n_image_handle in inputs["data"]:
            n_image = n_image_handle.get()
            data_id = n_image_handle.dataId
            band = str(data_id.band.name)
            patch = int(data_id.patch.id)
            median = np.nanmedian(n_image.array)
            stdev = np.nanstd(n_image.array)

            bands.append(band)
            patches.append(patch)
            medians.append(median)
            stdevs.append(stdev)

            band_patch_stats = []
            for threshold in self.config.threshold_list:
                # Calculate the percentage of the image with an image depth
                # above the given threshold.
                stat = np.sum(n_image.array > threshold) * 100 / (n_image.getHeight() * n_image.getWidth())
                band_patch_stats.append(stat)

            stats.append(band_patch_stats)

            # Calculate the quantiles for image depth
            # across the whole n_image array.
            quantile = list(np.percentile(n_image.array, q=self.config.quantile_list))
            quantiles.append(quantile)

        threshold_col_names = [
            f"depth_above_threshold_{threshold}" for threshold in self.config.threshold_list
        ]
        quantile_col_names = [f"depth_{q}_percentile" for q in self.config.quantile_list]

        # Construct the Astropy table
        data = [patches, bands, medians, stdevs] + list(zip(*stats)) + list(zip(*quantiles))
        names = ["patch", "band", "medians", "stdevs"] + threshold_col_names + quantile_col_names
        dtype = (
            ["int", "str", "float", "float"]
            + ["float" for x in range(len(list(zip(*stats))))]
            + ["int" for y in range(len(list(zip(*quantiles))))]
        )
        t = Table(data=data, names=names, dtype=dtype)
        return Struct(statTable=t)
