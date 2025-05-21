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
    "NImageCoaddSummaryAnalysisConfig",
    "NImageCoaddSummaryAnalysisTask",
)
from lsst.pex.config import ListField
from lsst.pipe.base import InputQuantizedConnection, OutputQuantizedConnection, QuantumContext
from lsst.pipe.base import connectionTypes as cT
from lsst.pipe.base import Struct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class NImageCoaddSummaryAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("tract", "skymap"),
    defaultTemplates={"inputName": "template_coadd_n_image", "outputName": "template_coadd_n_image_table"},
):
    data = cT.Input(
        doc="Number of input images per pixel summary statistics to load from the butler.",
        name="{inputName}",
        storageClass="ImageU",
        multiple=True,
        dimensions=("tract", "patch", "band", "skymap"),
        deferLoad=True,
    )

    statTable = cT.Output(
        doc="Table with n_image stats.",
        name="{outputName}",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )


class NImageCoaddSummaryAnalysisConfig(AnalysisBaseConfig, pipelineConnections=NImageCoaddSummaryAnalysisConnections):
    threshold_list = ListField(
        default=[1, 3, 5, 12],
        dtype=int,
        doc="The n_image pixel value thresholds.",
    )


class NImageCoaddSummaryAnalysisTask(AnalysisPipelineTask):
    ConfigClass = NImageCoaddSummaryAnalysisConfig
    _DefaultName = "nImageCoaddSummaryAnalysis"

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
                stat = np.sum(n_image.array > threshold)/(n_image.getHeight()*n_image.getWidth())*100
                band_patch_stats.append(stat)

            stats.append(band_patch_stats)

        data = [patches, bands, medians, stdevs] + list(zip(*stats))
        names = ["patch", "band", "medians", "stdevs"] + [f"above_thresh_{threshold}" for threshold in self.config.threshold_list]
        t = Table(data=data, names=names)
        return Struct(statTable=t)
