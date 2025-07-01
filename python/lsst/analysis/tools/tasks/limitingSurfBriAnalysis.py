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
    "LimitingSurfBriTask",
    "LimitingSurfBriConfig",
    "LimitingSurfBriConnections",
    "LimitingSurfBriAnalysisConnections",
    "LimitingSurfBriAnalysisConfig",
    "LimitingSurfBriAnalysisTask",
)

import logging

import numpy as np

import lsst.afw.math as afwMath
from lsst.pex.config import Field, ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections
from lsst.pipe.base.connectionTypes import Input, Output

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

log = logging.getLogger(__name__)


class LimitingSurfBriConnections(PipelineTaskConnections, dimensions=("instrument", "visit")):
    """Connections class for LimitingSurfBriTask."""

    recalibrated_star_detectors = Input(
        name="recalibrated_star_detector",
        storageClass="ArrowAstropy",
        doc="Visit-based source table",
        multiple=True,
        dimensions=("instrument", "visit", "detector",),
        deferLoad=True,
    )

    photoCalib = Input(
        name="visit_image.photoCalib",
        storageClass="PhotoCalib",
        doc="Photometric calibration associated with the visit image.",
        multiple=True,
        dimensions=("instrument", "visit", "detector",),
        deferLoad=True,
    )

    wcs = Input(
        name="visit_image.wcs",
        storageClass="wcs",
        doc="WCS header associated with the visit image.",
        multiple=True,
        dimensions=("instrument", "visit", "detector",),
        deferLoad=True,
    )

    limiting_surface_brightness_hist = Output(
        name="limiting_surface_brightness_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the "
        "aggregated limiting surface brightness dataset, i.e. limiting surface brightnesses estimated for "
        " all detectors in a single visit.",
        dimensions=("instrument", "visit")
    )


class LimitingSurfBriConfig(PipelineTaskConfig,
                            pipelineConnections=LimitingSurfBriConnections):
    """Config class for LimitingSurfBriTask."""

    bin_range = ListField[float](
        doc="The lower and upper range for the histogram bins, in ABmag.",
        default=[32, 26],
    )
    bin_width = Field[float](
        doc="The width of each histogram bin, in ABmag.",
        default=0.1,
    )


class LimitingSurfBriTask(PipelineTask):
    """A task for measuring the 3sigma limiting surface brightness on 10
    arcsecond scales for a given visit image.  This is the current standard
    depth metric within the low-surface-brightness science community.

    Reference: Roman, J., Trujillo, I., & Montes, M., 2020, A & A, 644, 42
    """

    ConfigClass = LimitingSurfBriConfig
    _DefaultName = "limitingSurfBri"

    def __init__(self, initInputs=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        limiting_surface_brightness_hist = self.run(**{k: v for k, v in inputs.items()})
        butlerQC.put(limiting_surface_brightness_hist, outputRefs.delta_skyCorr_hist)

    def run(self, recalibrated_star_detector, photoCalib, wcs):
        # Generate lookup tables for per-detector catalogues
        lookup_recalibrated_star_detector = {x.dataId: x for x in recalibrated_star_detector}
        lookup_photoCalib = {x.dataId: x for x in photoCalib}
        lookup_wcs = {x.dataId: x for x in wcs}

        # Set up the global histogram.
        bin_edges = np.arange(
            self.config.bin_range[0],
            self.config.bin_range[1] + self.config.bin_width,
            self.config.bin_width,
        )
        hist = np.zeros(len(bin_edges) - 1)
        log.info("Generating a histogram containing %d bins.", len(hist))

        # Loop over the recalibrated_star_detector data
        for dataId in lookup_recalibrated_star_detector.keys():
            # Get the detector catalogue
            recalibrated_star_catalogue = lookup_recalibrated_star_detector[dataId].get()
            # And the photometric calibration
            nanojanskyToInstFlux = lookup_photoCalib[dataId].get().nanojanskyToInstFlux(1)
            instFluxToMagnitude = lookup_photoCalib[dataId].get().instFluxToMagnitude(1)
            # And the pixel scale
            pxScale = lookup_wcs[dataId].get().getPixelScale().asArcseconds()

            # Isolate the sky sources, 9px radius apertures only
            isSky = (recalibrated_star_catalogue["sky_source"] > 0)
            skySources = recalibrated_star_catalogue[isSky]["ap09Flux"]

            # Derive the clipped standard deviation of sky sources in nJy
            ctrl = afwMath.StatisticsControl(3, 3)  # TODO: make these configurable
            ctrl.setNanSafe(True)
            statistic = afwMath.stringToStatisticsProperty("STDEVCLIP")
            sigSkySources = afwMath.makeStatistics(skySources, statistic, ctrl).getValue()

            # Derive limiting surface brightness.  3sigma, on 10"x10" scales
            nPix = np.pi*9**2  # Number of pixels within the circular aperture
            sigSkySources *= nanojanskyToInstFlux
            pixScaleRatio = np.sqrt(pxScale**2 / (nPix*pxScale**2))
            sigma = sigSkySources / pixScaleRatio
            muLim = -2.5*np.log10((3*sigma)/(pxScale*10)) + instFluxToMagnitude

            # Compute the per-detector histogram; update the global histogram.
            hist_det, _ = np.histogram(muLim, bins=bin_edges)
            hist += hist_det

        # Return results.
        num_populated_bins = len([x for x in hist if x == 0])
        log.info("Populated %d of %d histogram bins.", len(hist) - num_populated_bins, len(hist))
        bin_mid = bin_edges[:-1] + (self.config.bin_width / 2)
        limiting_surface_brightness_hist = dict(
            hist=hist, bin_lower=bin_edges[:-1], bin_upper=bin_edges[1:], bin_mid=bin_mid
        )
        return limiting_surface_brightness_hist


class LimitingSurfBriAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument", "visit"),
    defaultTemplates={"outputName": "limitingSurfBri"},
):
    data = Input(
        name="limiting_surface_brightness_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the "
        "aggregated limiting surface brightness dataset, i.e. limiting surface brightnesses estimated for "
        " all detectors in a single visit.",
        deferLoad=True,
        dimensions=("instrument", "visit"),
    )


class LimitingSurfBriAnalysisConfig(AnalysisBaseConfig,
                                    pipelineConnections=LimitingSurfBriAnalysisConnections):
    pass


class LimitingSurfBriAnalysisTask(AnalysisPipelineTask):
    ConfigClass = LimitingSurfBriAnalysisConfig
    _DefaultName = "limitingSurfBriAnalysis"
