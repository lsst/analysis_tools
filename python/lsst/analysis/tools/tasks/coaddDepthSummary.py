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
from astropy.table import Table, join

from lsst.pex.config import ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base import connectionTypes as cT


class CoaddDepthSummaryConnections(
    PipelineTaskConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"coaddName": ""},  # set as either deep or template in the pipeline
):
    mask_data = cT.Input(
        doc="Coadd to load from the butler.",
        name="{coaddName}_coadd.mask",
        storageClass="MaskX",
        multiple=True,
        dimensions=("tract", "patch", "band", "skymap"),
        deferLoad=True,
    )

    n_image_data = cT.Input(
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
        default=[1, 2, 3, 5, 12],
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
        # Calculate coadd metrics.
        coadd_bands = []
        coadd_patches = []
        chip_gap_percents = []

        for mask_handle in inputs["mask_data"]:
            mask = mask_handle.get()
            data_id = mask_handle.dataId
            band = str(data_id.band.name)
            patch = int(data_id.patch.id)

            coadd_bands.append(band)
            coadd_patches.append(patch)

            no_data = mask.getPlaneBitMask("NO_DATA")
            rejected = mask.getPlaneBitMask("REJECTED")
            mask_array = mask.array

            # Calculate pixels that are NO_DATA but not REJECTED.
            no_data_flag = (mask_array & no_data) != 0
            rejected_flag = (mask_array & rejected) != 0

            chip_gap_percent = ((no_data_flag & ~rejected_flag).sum() * 100
                                / (mask.getHeight() * mask.getWidth()))

            chip_gap_percents.append(chip_gap_percent)

        # Construct the Astropy table for coadd information.
        data = [coadd_patches, coadd_bands, chip_gap_percents]
        names = ["patch", "band", "chip_gap_percent"]
        dtype = ["int", "str", "float"]
        coadd_table = Table(data=data, names=names, dtype=dtype)

        # Calculate n_image metrics.
        n_image_bands = []
        n_image_patches = []
        means = []
        medians = []
        stdevs = []
        stats = []
        quantiles = []

        for n_image_handle in inputs["n_image_data"]:
            n_image = n_image_handle.get()
            data_id = n_image_handle.dataId
            band = str(data_id.band.name)
            patch = int(data_id.patch.id)
            mean = np.nanmean(n_image.array)
            median = np.nanmedian(n_image.array)
            stdev = np.nanstd(n_image.array)

            n_image_bands.append(band)
            n_image_patches.append(patch)
            means.append(mean)
            medians.append(median)
            stdevs.append(stdev)

            band_patch_stats = []
            for threshold in self.config.threshold_list:
                # Calculate the percentage of the image with an image depth
                # above the given threshold.
                stat = np.sum(n_image.array >= threshold) * 100 / (n_image.getHeight() * n_image.getWidth())
                band_patch_stats.append(stat)

            stats.append(band_patch_stats)

            # Calculate the quantiles for image depth across n_image array.
            quantile = list(np.percentile(n_image.array, q=self.config.quantile_list))
            quantiles.append(quantile)

        threshold_col_names = [
            f"depth_above_threshold_{threshold}" for threshold in self.config.threshold_list
        ]
        quantile_col_names = [f"depth_{q}_percentile" for q in self.config.quantile_list]

        # Construct the Astropy table for n_image information.
        data = (
            [n_image_patches, n_image_bands, means, medians, stdevs]
            + list(zip(*stats))
            + list(zip(*quantiles))
        )
        names = (
            ["patch", "band", "mean", "median", "stdevs"]
            + threshold_col_names
            + quantile_col_names
        )
        dtype = (
            ["int", "str", "float", "float", "float"]
            + ["float" for x in range(len(list(zip(*stats))))]
            + ["int" for y in range(len(list(zip(*quantiles))))]
        )
        n_image_table = Table(data=data, names=names, dtype=dtype)

        # Combine tables.
        combined_table = join(coadd_table, n_image_table, keys=["patch", "band"])
        return Struct(statTable=combined_table)
