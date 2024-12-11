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

from lsst.pipe.base import InputQuantizedConnection, OutputQuantizedConnection, QuantumContext
from lsst.pipe.base import connectionTypes as cT
from lsst.pipe.base import Struct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class NImageCoaddSummaryAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("tract", "band", "skymap"),
    defaultTemplates={"inputName": "goodSeeingCoadd_nImage", "outputName": "nImageTable"},
):
    data = cT.Input(
        doc="Number of input images per pixel summary statistics to load from the butler",
        name="{inputName}",
        storageClass="ImageU",
        multiple=True,
        dimensions=("tract", "patch", "band", "skymap"),
        deferLoad=True,
    )

    statTable = cT.Output(
        doc="Table with n_image stats",
        name="{outputName}_tract",
        storageClass="ArrowAstropy",
        dimensions=("tract", "skymap"),
    )


class NImageCoaddSummaryAnalysisConfig(AnalysisBaseConfig, pipelineConnections=NImageCoaddSummaryAnalysisConnections):
    pass


class NImageCoaddSummaryAnalysisTask(AnalysisPipelineTask):
    ConfigClass = NImageCoaddSummaryAnalysisConfig
    _DefaultName = "nImageCoaddSummaryAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(inputs)
        butlerQC.put(outputs, outputRefs)
    
    def run(self, inputs):
        bands = []
        patches = []
        stats = []
        for n_image_handle in inputs["data"]:
            data_id = n_image_handle.dataId
            band = str(data_id.band.name)
            patch = int(data_id.patch.id)

            bands.append(band)
            patches.append(patch)
            n_image = n_image_handle.get()
            stat = np.sum(n_image.array > 3)/(n_image.getHeight()*n_image.getWidth())*100
            stats.append(stat)

        t = Table(data=[patches, bands, stats], names=["Patch", "Band", "Stat"])
        
        return Struct(statTable=t)
        
        # inputs["num_initial_bgs"] = len(inputs["calexpBackgrounds"][0].get())
        # delta_skyCorr_hist = self.run(**{k: v for k, v in inputs.items() if k != "calexpBackgrounds"})
        # butlerQC.put(delta_skyCorr_hist, outputRefs.delta_skyCorr_hist)

    # def runQuantum(
    #     self,
    #     butlerQC: QuantumContext,
    #     inputRefs: InputQuantizedConnection,
    #     outputRefs: OutputQuantizedConnection,
    # ) -> None:
    #     # Docstring inherited.

    #     inputs = butlerQC.get(inputRefs)

    #     summary = inputs["data"].__dict__

    #     outputs = self.run(data=summary)
    #     butlerQC.put(outputs, outputRefs)
