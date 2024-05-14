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

__all__ = (
    "AllDiaSourceHistTask",
    "AllDiaSourceHistConfig",
    "AllDiaSourceHistConnections",
    "AllDiaSourceAnalysisConnections",
    "AllDiaSourceAnalysisConfig",
    "AllDiaSourceAnalysisTask",
)

import logging

import numpy as np
from lsst.pex.config import Field, ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections
from lsst.pipe.base.connectionTypes import Input, Output

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

log = logging.getLogger(__name__)


class AllDiaSourceHistConnections(PipelineTaskConnections, dimensions=("instrument", "skymap"), defaultTemplates={"coaddName": "goodSeeing", "fakesType": ""}):
    """Connections class for AllDiaSourceHistTask."""
    
    sourceTables = Input(
                 doc="CcdVisit-based DiaSource table to load from the butler",
                 name="{fakesType}{coaddName}Diff_assocDiaSrc",
                 storageClass="DataFrame",
                 deferLoad=True,
                 multiple=True,
                 dimensions=("visit", "band", "detector"),
    )
    
    flux_hist = Output(
        name="flux_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the "
        "aggregated DIA source flux.",
        dimensions=("instrument", "skymap"),
    )


class AllDiaSourceHistConfig(PipelineTaskConfig, pipelineConnections=AllDiaSourceHistConnections):
    """Config class for AllDiaSourceHistTask."""

    bin_range = ListField[float](
        doc="The lower and upper range for the histogram bins, in nJy.",
        default=[-3e3, 3e3],
    )
    flux_type = Field[str](
        doc="The type of flux (column name).",
        default='psfFlux',
    )
    num_bins = Field[int](
        doc="The number of bins.",
        default=200,
    )


class AllDiaSourceHistTask(PipelineTask):
    """A task for generating a histogram of counts in the specified flux of DIA Sources.
    """

    ConfigClass = AllDiaSourceHistConfig
    _DefaultName = "allDiaSourceHist"

    def __init__(self, initInputs=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        flux_hist = self.run(inputs["sourceTables"])
        breakpoint()
        butlerQC.put(flux_hist, outputRefs.flux_hist)

        # inputs = butlerQC.get(inputRefs)
        # inputs["num_initial_bgs"] = len(inputs["calexpBackgrounds"][0].get())
        # delta_skyCorr_hist = self.run(**{k: v for k, v in inputs.items() if k != "calexpBackgrounds"})
        # butlerQC.put(delta_skyCorr_hist, outputRefs.delta_skyCorr_hist)

    def run(self, sourceTables):
        """Generate a histogram of counts in the ...

        Parameters
        ----------
        sourceTables : `list`[`~lsst.daf.butler.DeferredDatasetHandle`]
            CcdVisit-based DiaSource table to load from the butler.

        Returns
        -------
        flux_hist : `dict`[`str`, `~numpy.ndarray`]
            A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the 
            aggregated DIA source flux.
        """
        # Generate lookup tables for the DIA source tables.
        lookup_sourceTables = {x.dataId: x for x in sourceTables}

        # Set up the global histogram.
        bin_edges = np.linspace(
            self.config.bin_range[0],
            self.config.bin_range[1],
            self.config.num_bins,
        )
        hist = np.zeros(self.config.num_bins)
        log.info("Generating a histogram containing %d bins.", len(hist))

        # Loop over the skyCorr/injected_skyCorr data.
        for dataId in lookup_sourceTables.keys():
            # Get the skyCorr/injected_skyCorr data.
            sourceTable = lookup_sourceTables[dataId].get()

            # Compute the per-detector histogram; update the global histogram.
            finiteData = np.isfinite(sourceTable[self.config.flux_type])
            hist_det, _ = np.histogram(sourceTable[self.config.flux_type][finiteData], bins=self.config.num_bins)
            hist += hist_det

        # Return results.
        num_populated_bins = len([x for x in hist if x == 0])
        log.info("Populated %d of %d histogram bins.", len(hist) - num_populated_bins, len(hist))
        # bin_mid = bin_edges[:-1] + (self.config.bin_width / 2)
        bin_mid = bin_edges[int(len(bin_edges) / 2)]
        flux_hist = dict(
            hist=hist-1, bin_lower=bin_edges, bin_upper=bin_edges, bin_mid=bin_mid
        )
        return flux_hist


class AllDiaSourceAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "skymap"),
    defaultTemplates={"outputName": "AllDiaSource"},
):
    data = Input(
        name="flux_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges "
        "for the selected flux type for all DIA sources.",
        deferLoad=True,
        dimensions=("instrument", "skymap"),
    )


class AllDiaSourceAnalysisConfig(AnalysisBaseConfig, pipelineConnections=AllDiaSourceAnalysisConnections):
    pass


class AllDiaSourceAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AllDiaSourceAnalysisConfig
    _DefaultName = "allDiaSourceAnalysis"
