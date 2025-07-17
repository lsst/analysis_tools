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

import lsst.afw.math as afwMath
import numpy as np
from astropy.table import Table
from lsst.pex.config import Field, ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base.connectionTypes import Input, Output
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

log = logging.getLogger(__name__)


class LimitingSurfBriConnections(
    PipelineTaskConnections,
    dimensions=(),
    defaultTemplates={
        "detectionTableName": "",
        "photoCalibName": "",
        "wcsName": "",
    },
):
    """Connections class for LimitingSurfBriTask."""

    data = Input(
        doc="Visit- or coadd-level object table",
        name="{detectionTableName}",
        storageClass="ArrowAstropy",
        multiple=True,
        dimensions=(),
        deferLoad=True,
    )

    skymap = Input(
        doc="The skymap that covers the originating data's tract.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    photoCalib = Input(
        name="{photoCalibName}",
        storageClass="PhotoCalib",
        doc="Photometric calibration associated with the originating data image.",
        multiple=True,
        dimensions=(),
        deferLoad=True,
    )

    wcs = Input(
        name="{wcsName}",
        storageClass="Wcs",
        doc="WCS header associated with the originating data image.",
        multiple=True,
        dimensions=(),
        deferLoad=True,
    )

    limiting_surface_brightness_hist = Output(
        name="limiting_surface_brightness_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the "
        "aggregated limiting surface brightness dataset, i.e. limiting surface brightnesses estimated for "
        " all detectors in a single visit or all patches in a single tract.",
        dimensions=(),
    )

    limiting_surface_brightness_table = Output(
        name="limiting_surface_brightness_table",
        storageClass="ArrowAstropy",
        doc="A table containing two columns: the detector or patch IDs and the values of limiting surface "
        "brightness derived for those detectors or patches.",
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        self.data = Input(
            doc=self.data.doc,
            name=self.data.name,
            storageClass=self.data.storageClass,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputTableDimensions)),
            multiple=self.data.multiple,
        )
        self.photoCalib = Input(
            doc=self.photoCalib.doc,
            name=self.photoCalib.name,
            storageClass=self.photoCalib.storageClass,
            deferLoad=self.photoCalib.deferLoad,
            dimensions=frozenset(sorted(config.inputCalibDimensions)),
            multiple=self.photoCalib.multiple,
        )
        self.wcs = Input(
            doc=self.wcs.doc,
            name=self.wcs.name,
            storageClass=self.wcs.storageClass,
            deferLoad=self.wcs.deferLoad,
            dimensions=frozenset(sorted(config.inputCalibDimensions)),
            multiple=self.wcs.multiple,
        )
        self.limiting_surface_brightness_hist = Output(
            doc=self.limiting_surface_brightness_hist.doc,
            name=self.limiting_surface_brightness_hist.name,
            storageClass=self.limiting_surface_brightness_hist.storageClass,
            dimensions=frozenset(sorted(config.outputDataDimensions)),
        )

        self.limiting_surface_brightness_table = Output(
            doc=self.limiting_surface_brightness_table.doc,
            name=self.limiting_surface_brightness_table.name,
            storageClass=self.limiting_surface_brightness_table.storageClass,
            dimensions=frozenset(sorted(config.outputDataDimensions)),
        )

        assert config is not None, "Missing required config object."

        if "tract" not in config.inputTableDimensions:
            del self.skymap

        self.dimensions.update(frozenset(sorted(config.outputDataDimensions)))


class LimitingSurfBriConfig(PipelineTaskConfig, pipelineConnections=LimitingSurfBriConnections):
    """Config class for LimitingSurfBriTask."""

    inputTableDimensions = ListField[str](
        doc="Dimensions of the input object table data.",
        default=(),
        optional=False,
    )
    inputCalibDimensions = ListField[str](
        doc="Dimensions of the input calibration data.",
        default=(),
        optional=False,
    )
    outputDataDimensions = ListField[str](
        doc="Dimensions of the output histogram and table data.",
        default=(),
        optional=False,
    )

    bin_range = ListField[float](
        doc="The lower and upper range for the histogram bins, in ABmag.",
        default=[26.5, 32],
    )
    bin_width = Field[float](
        doc="The width of each histogram bin, in ABmag.",
        default=0.1,
    )


class LimitingSurfBriTask(PipelineTask):
    """A task for measuring the 3sigma limiting surface brightness on 10
    arcsecond scales for a given image.  This is currently a widely accepted
    metric for depth within the low surface brightness community.

    Reference: Roman, J., Trujillo, I., & Montes, M., 2020, A & A, 644, 42
    """

    ConfigClass = LimitingSurfBriConfig
    _DefaultName = "limitingSurfBri"

    def __init__(self, initInputs=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        limiting_surface_brightness_struct = self.run(**{k: v for k, v in inputs.items() if k != "skymap"})
        butlerQC.put(limiting_surface_brightness_struct, outputRefs)

    def run(self, data, photoCalib, wcs):
        """Output a histogram and table of per-detector or per-patch limiting
        surface brightnesses for the input visit or tract coadd

        Parameters
        ----------
        data : `list`
            List of dicts with information respecting the extracted source
            detection catalogues
        photoCalib : `list`
            List of dicts with information respecting the extracted image
            photometric calibrations
        wcs : `list`
            List of dicts with information respecting the extracted image
            world coordinate systems

        Returns
        -------
        `pipe.base.Struct` containing `dict` and `astropy.table.Table`
            Dictionary containing limiting surface brightness histogram
            bin parameters and counts; table containing per-detector or per-
            patch limiting surface brightness values
        """
        # Generate lookup tables
        lookup_photoCalib = {x.dataId: x for x in photoCalib}
        lookup_wcs = {x.dataId: x for x in wcs}

        # Retrieve the source catalogue and photometric band if required
        source_catalogue = data[0].get()
        if "band" in self.config.inputCalibDimensions:
            band = photoCalib[0].dataId["band"] + "_"
        else:
            band = ""

        # Set up the global histogram.
        bin_edges = np.arange(
            self.config.bin_range[0],
            self.config.bin_range[1] + self.config.bin_width,
            self.config.bin_width,
        )
        hist = np.zeros(len(bin_edges) - 1)
        log.info("Generating a histogram containing %d bins.", len(hist))

        # Set up the global table
        limiting_surface_brightness_table = Table(
            names=["imageID", "limiting_surface_brightness"],
            dtype=[np.int64, np.float64],
        )

        # Loop over the detector or patch calibration data
        for dataId in lookup_photoCalib.keys():
            # And the photometric calibration
            instFluxToMagnitude = lookup_photoCalib[dataId].get().instFluxToMagnitude(1)
            # And the pixel scale
            pxScale = lookup_wcs[dataId].get().getPixelScale().asArcseconds()

            # Isolate the sky sources per image, 9px radius apertures
            idKey = "detector" if "detector" in dataId else "patch"
            skyKey = "sky_source" if "detector" in dataId else "sky_object"
            isImage = source_catalogue[idKey] == dataId[idKey]
            isSky = source_catalogue[skyKey] > 0
            skySources = source_catalogue[isImage & isSky][band + "ap09Flux"]

            # Derive the clipped standard deviation of sky sources in nJy
            nPix = np.pi * 9**2  # Number of pixels within the circular aperture
            ctrl = afwMath.StatisticsControl(3, 3)
            ctrl.setNanSafe(True)
            statistic = afwMath.stringToStatisticsProperty("STDEVCLIP")
            sigSkySources = afwMath.makeStatistics(skySources / nPix, statistic, ctrl).getValue()

            # Derive limiting surface brightness.  3sigma, on 10"x10" scales
            pixScaleRatio = np.sqrt(pxScale**2 / (nPix * pxScale**2))
            sigma = sigSkySources / pixScaleRatio
            muLim = -2.5 * np.log10((3 * sigma) / (pxScale * 10)) + instFluxToMagnitude

            # Compute the histogram; update the global histogram.
            hist_det, _ = np.histogram(muLim, bins=bin_edges)
            hist += hist_det

            # Append a new row to the table
            limiting_surface_brightness_table.add_row([dataId[idKey], muLim])

        # Return results.
        num_populated_bins = len([x for x in hist if x == 0])
        log.info("Populated %d of %d histogram bins.", len(hist) - num_populated_bins, len(hist))
        bin_mid = bin_edges[:-1] + (self.config.bin_width / 2)
        limiting_surface_brightness_hist = dict(
            hist=hist, bin_lower=bin_edges[:-1], bin_upper=bin_edges[1:], bin_mid=bin_mid
        )
        return Struct(
            limiting_surface_brightness_hist=limiting_surface_brightness_hist,
            limiting_surface_brightness_table=limiting_surface_brightness_table,
        )


class LimitingSurfBriAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=(),
    defaultTemplates={"outputName": "limitingSurfBri"},
):
    data = Input(
        name="limiting_surface_brightness_hist",
        storageClass="ArrowNumpyDict",
        doc="A dictionary containing the histogram values, bin mid points, and bin lower/upper edges for the "
        "aggregated limiting surface brightness dataset, i.e. limiting surface brightnesses estimated for "
        " all detectors in a single visit or all patches in a single tract.",
        deferLoad=True,
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        self.data = Input(
            name=self.data.name,
            storageClass=self.data.storageClass,
            doc=self.data.doc,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputDataDimensions)),
        )


class LimitingSurfBriAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=LimitingSurfBriAnalysisConnections
):
    inputDataDimensions = ListField[str](
        doc="Dimensions of the input histogram and table data.",
        default=(),
        optional=False,
    )


class LimitingSurfBriAnalysisTask(AnalysisPipelineTask):
    ConfigClass = LimitingSurfBriAnalysisConfig
    _DefaultName = "limitingSurfBriAnalysis"
