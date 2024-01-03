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
import logging

import numpy as np
from lsst.pex.config import Field, ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections
from lsst.pipe.base.connectionTypes import Input, Output

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

log = logging.getLogger(__name__)


class DeltaSkyCorrHistConnections(PipelineTaskConnections, dimensions=("instrument", "visit")):
    """Connections class for DeltaSkyCorrHistTask."""

    skyCorrs = Input(
        name="skyCorr",
        storageClass="Background",
        doc="Sky correction background models from a run without any synthetic source injection.",
        multiple=True,
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
    )
    injected_skyCorrs = Input(
        name="injected_skyCorr",
        storageClass="Background",
        doc="Sky correction background models from a run with synthetic sources injected into the data.",
        multiple=True,
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
    )
    calexpBackgrounds = Input(
        name="calexpBackground",
        storageClass="Background",
        doc="Initial per-detector background models associated with the calibrated exposure.",
        multiple=True,
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
    )
    photoCalib = Input(
        name="calexp.photoCalib",
        storageClass="PhotoCalib",
        doc="Photometric calibration associated with the calibrated exposure.",
        multiple=True,
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
    )
    delta_skyCorr_hist = Output(
        name="delta_skyCorr_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the "
        "aggregated skyCorr difference dataset, i.e., the difference between the injected and non-injected "
        "sky correction background models.",
        dimensions=("instrument", "visit"),
    )


class DeltaSkyCorrHistConfig(PipelineTaskConfig, pipelineConnections=DeltaSkyCorrHistConnections):
    """Config class for DeltaSkyCorrHistTask."""

    bin_range = ListField[float](
        doc="The lower and upper range for the histogram bins, in nJy.",
        default=[-1, 1],
    )
    bin_width = Field[float](
        doc="The width of each histogram bin, in nJy.",
        default=0.0001,
    )


class DeltaSkyCorrHistTask(PipelineTask):
    """A task for generating a histogram of counts in the difference image
    between an injected sky correction frame and a non-injected sky correction
    frame (i.e., injected_skyCorr - skyCorr).
    """

    ConfigClass = DeltaSkyCorrHistConfig
    _DefaultName = "deltaSkyCorrHist"

    def __init__(self, initInputs=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["num_initial_bgs"] = len(inputs["calexpBackgrounds"][0].get())
        delta_skyCorr_hist = self.run(**{k: v for k, v in inputs.items() if k != "calexpBackgrounds"})
        butlerQC.put(delta_skyCorr_hist, outputRefs.delta_skyCorr_hist)

    def run(self, skyCorrs, injected_skyCorrs, num_initial_bgs, photoCalib):
        """Generate a histogram of counts in the difference image between an
        injected sky correction frame and a non-injected sky correction frame
        (i.e., injected_skyCorr - skyCorr).

        Parameters
        ----------
        skyCorrs : `list`[`~lsst.daf.butler.DeferredDatasetHandle`]
            Sky correction background models from a run without any synthetic
            source injection.
            These deferred dataset handles should normally resolve to
            `~lsst.afw.math.BackgroundList` objects.
        injected_skyCorrs : `list`[`~lsst.daf.butler.DeferredDatasetHandle`]
            Sky correction background models from a run with synthetic sources
            injected into the data.
            These deferred dataset handles should normally resolve to
            `~lsst.afw.math.BackgroundList` objects.
        num_initial_bgs : `int`
            The length of the initial per-detector background model list.
            This number of background models will be skipped from the start of
            each skyCorr/injected_skyCorr background model list.
            See the Notes section for more details.
        photoCalib : `list`[`~lsst.daf.butler.DeferredDatasetHandle`]
            Photometric calibration, for conversion from counts to nJy.

        Returns
        -------
        delta_skyCorr_hist : `dict`[`str`, `~numpy.ndarray`]
            A dictionary containing the histogram values and bin lower/upper
            edges for the skyCorr difference dataset.

        Notes
        -----
        The first N background elements in the skyCorr/injected_skyCorr
        background list are the inverse of the initial per-detector background
        solution.
        The effect of this is that adding a sky correction frame to a
        background-subtracted calibrated exposure will undo the per-detector
        background solution and apply the full focal plane sky correction in
        its place.

        For this task, we only want to compare the extra (subtractive) sky
        correction components, so we skip the first N background models from
        the sky frame.
        """
        # Generate lookup tables for the skyCorr/injected_skyCorr data.
        lookup_skyCorrs = {x.dataId: x for x in skyCorrs}
        lookup_injected_skyCorrs = {x.dataId: x for x in injected_skyCorrs}
        lookup_photoCalib = {x.dataId: x for x in photoCalib}

        # Set up the global histogram.
        bin_edges = np.arange(
            self.config.bin_range[0],
            self.config.bin_range[1] + self.config.bin_width,
            self.config.bin_width,
        )
        hist = np.zeros(len(bin_edges) - 1)
        log.info("Generating a histogram containing %d bins.", len(hist))

        # Loop over the skyCorr/injected_skyCorr data.
        for dataId in lookup_injected_skyCorrs.keys():
            # Get the skyCorr/injected_skyCorr data.
            skyCorr = lookup_skyCorrs[dataId].get()
            injected_skyCorr = lookup_injected_skyCorrs[dataId].get()
            # And the photometric calibration
            instFluxToNanojansky = lookup_photoCalib[dataId].get().instFluxToNanojansky(1)

            # Isolate the extra (subtractive) sky correction components.
            skyCorr_extras = skyCorr.clone()
            skyCorr_extras._backgrounds = skyCorr_extras._backgrounds[num_initial_bgs:]
            injected_skyCorr_extras = injected_skyCorr.clone()
            injected_skyCorr_extras._backgrounds = injected_skyCorr_extras._backgrounds[num_initial_bgs:]

            # Create the delta_skyCorr array.
            delta_skyCorr_det = injected_skyCorr_extras.getImage().array - skyCorr_extras.getImage().array
            delta_skyCorr_det *= instFluxToNanojansky  # Convert image to nJy

            # Compute the per-detector histogram; update the global histogram.
            hist_det, _ = np.histogram(delta_skyCorr_det, bins=bin_edges)
            hist += hist_det

        # Return results.
        num_populated_bins = len([x for x in hist if x == 0])
        log.info("Populated %d of %d histogram bins.", len(hist) - num_populated_bins, len(hist))
        bin_mid = bin_edges[:-1] + (self.config.bin_width / 2)
        delta_skyCorr_hist = dict(
            hist=hist, bin_lower=bin_edges[:-1], bin_upper=bin_edges[1:], bin_mid=bin_mid
        )
        return delta_skyCorr_hist


class DeltaSkyCorrAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "visit"),
    defaultTemplates={"outputName": "deltaSkyCorr"},
):
    data = Input(
        name="delta_skyCorr_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the "
        "aggregated skyCorr difference dataset, i.e., the difference between the injected and non-injected "
        "sky correction background models.",
        deferLoad=True,
        dimensions=("instrument", "visit"),
    )


class DeltaSkyCorrAnalysisConfig(AnalysisBaseConfig, pipelineConnections=DeltaSkyCorrAnalysisConnections):
    pass


class DeltaSkyCorrAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DeltaSkyCorrAnalysisConfig
    _DefaultName = "deltaSkyCorrAnalysis"
