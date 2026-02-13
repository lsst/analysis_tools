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
    "LimitingSurfaceBrightnessTask",
    "LimitingSurfaceBrightnessConfig",
    "LimitingSurfaceBrightnessConnections",
    "LimitingSurfaceBrightnessAnalysisConnections",
    "LimitingSurfaceBrightnessAnalysisConfig",
    "LimitingSurfaceBrightnessAnalysisTask",
)

import logging

import numpy as np
from astropy.table import Table

from lsst.afw.math import StatisticsControl, makeStatistics, stringToStatisticsProperty
from lsst.pex.config import Field, ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base.connectionTypes import Input, Output
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

log = logging.getLogger(__name__)


class LimitingSurfaceBrightnessConnections(
    PipelineTaskConnections,
    dimensions=(),
    defaultTemplates={
        "detectionTableName": "object_all",
        "photoCalibName": "deep_coadd.photoCalib",
        "wcsName": "deep_coadd.wcs",
    },
):
    """Connections class for LimitingSurfaceBrightnessTask."""

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

    limiting_surface_brightness_table = Output(
        name="limiting_surface_brightness_table",
        storageClass="ArrowAstropy",
        doc="A table containing two columns: the detector or patch IDs and the values of limiting surface "
        "brightness derived for those detectors or patches.",
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        # Update output table name for configurable dimensions
        dimen = "_visit" if "visit" in config.inputTableDimensions else "_tract"
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
        self.limiting_surface_brightness_table = Output(
            doc=self.limiting_surface_brightness_table.doc,
            name=self.limiting_surface_brightness_table.name + dimen,
            storageClass=self.limiting_surface_brightness_table.storageClass,
            dimensions=frozenset(sorted(config.outputDataDimensions)),
        )

        assert config is not None, "Missing required config object."

        if "tract" not in config.inputTableDimensions:
            del self.skymap

        self.dimensions.update(frozenset(sorted(config.outputDataDimensions)))


class LimitingSurfaceBrightnessConfig(
    PipelineTaskConfig, pipelineConnections=LimitingSurfaceBrightnessConnections
):
    """Config class for LimitingSurfaceBrightnessTask."""

    inputTableDimensions = ListField[str](
        doc="Dimensions of the input object table data.",
        default=("skymap", "tract"),
        optional=False,
    )
    inputCalibDimensions = ListField[str](
        doc="Dimensions of the input calibration data.",
        default=("tract", "band"),
        optional=False,
    )
    outputDataDimensions = ListField[str](
        doc="Dimensions of the output table data.",
        default=("tract", "band"),
        optional=False,
    )
    apertureSize = Field[int](
        doc="The size of the sky objects photometry aperture.",
        default=9,
    )


class LimitingSurfaceBrightnessTask(PipelineTask):
    """A task for measuring the 3 sigma limiting surface brightness on 10
    arcsecond scales for a given image.  This is currently a widely accepted
    metric for depth within the low surface brightness community.

    Reference: Roman, J., Trujillo, I., & Montes, M., 2020, A & A, 644, 42
    """

    ConfigClass = LimitingSurfaceBrightnessConfig
    _DefaultName = "limitingSurfaceBrightness"

    def __init__(self, initInputs=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        limiting_surface_brightness_struct = self.run(**{k: v for k, v in inputs.items() if k != "skymap"})
        butlerQC.put(limiting_surface_brightness_struct, outputRefs)

    def run(self, data, photoCalib, wcs):
        """Output a table of per-detector or per-patch limiting surface
        brightnesses for the input visit or tract coadd

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
        `pipe.base.Struct` containing `astropy.table.Table`
            Table containing per-detector or per-patch limiting surface
            brightness values
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
            skySources = source_catalogue[isImage & isSky][f"{band}ap{self.config.apertureSize:02d}Flux"]
            # Some patches contain no detections
            if len(skySources) == 0:
                muLim = np.nan
                limiting_surface_brightness_table.add_row([dataId[idKey], muLim])
                continue

            # Derive the clipped standard deviation of sky sources in nJy
            nPix = np.pi * self.config.apertureSize**2  # Number of pixels within the circular aperture
            ctrl = StatisticsControl(3, 3)
            ctrl.setNanSafe(True)
            statistic = stringToStatisticsProperty("STDEVCLIP")
            sigSkySources = makeStatistics(skySources / nPix, statistic, ctrl).getValue()

            # Derive limiting surface brightness.  3sigma, on 10"x10" scales
            pixScaleRatio = np.sqrt(pxScale**2 / (nPix * pxScale**2))
            sigma = sigSkySources / pixScaleRatio
            muLim = -2.5 * np.log10((3 * sigma) / (pxScale * 10)) + instFluxToMagnitude

            # Append a new row to the table
            limiting_surface_brightness_table.add_row([dataId[idKey], muLim])

        # Return results.
        return Struct(limiting_surface_brightness_table=limiting_surface_brightness_table)


class LimitingSurfaceBrightnessAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=(),
    defaultTemplates={"outputName": "limitingSurfaceBrightness"},
):
    data = Input(
        name="limiting_surface_brightness_table",
        storageClass="ArrowAstropy",
        doc="A table containing two columns: the detector or patch IDs and the values of limiting surface "
        "brightness derived for those detectors or patches.",
        deferLoad=True,
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        # Update input table name for configurable dimensions
        dimen = "_visit" if "visit" in config.inputDataDimensions else "_tract"
        self.data = Input(
            name=self.data.name + dimen,
            storageClass=self.data.storageClass,
            doc=self.data.doc,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputDataDimensions)),
        )
        self.dimensions.update(frozenset(sorted(config.inputDataDimensions)))


class LimitingSurfaceBrightnessAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=LimitingSurfaceBrightnessAnalysisConnections
):
    inputDataDimensions = ListField[str](
        doc="Dimensions of the input table data.",
        default=(),
        optional=False,
    )


class LimitingSurfaceBrightnessAnalysisTask(AnalysisPipelineTask):
    ConfigClass = LimitingSurfaceBrightnessAnalysisConfig
    _DefaultName = "limitingSurfaceBrightnessAnalysis"
