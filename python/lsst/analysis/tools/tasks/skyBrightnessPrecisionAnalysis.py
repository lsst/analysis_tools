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
    "SkyBrightnessPrecisionTask",
    "SkyBrightnessPrecisionConfig",
    "SkyBrightnessPrecisionConnections",
    "SkyBrightnessPrecisionAnalysisConnections",
    "SkyBrightnessPrecisionAnalysisConfig",
    "SkyBrightnessPrecisionAnalysisTask",
)

import logging

import lsst.afw.math as afwMath
import numpy as np
from astropy.table import Table
from lsst.pex.config import Field, ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base.connectionTypes import Input, Output
from lsst.skymap import BaseSkyMap
from lsst.afw.geom import SpanSet


from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

log = logging.getLogger(__name__)


class SkyBrightnessPrecisionConnections(
    PipelineTaskConnections,
    defaultTemplates={
        "detectionTableName": "sourceTable",
        "calexpName": "calexp",
        "backgroundName": "calexpBackground",
        "photoCalibName": "calexp.photoCalib",
        "wcsName": "calexp.wcs",
    },
    dimensions=(),
):
    """Connections for the SkyBrightnessPrecisionTask."""

    objTable = Input(
        doc="Visit- or coadd-level object table",
        name="{detectionTableName}",
        storageClass="ArrowAstropy",
        multiple=True,
        dimensions=(),
        deferLoad=True,
    )

    calexp = Input(
        name="{calexpName}",
        storageClass="ExposureF",
        doc="Calibrated exposure (used for photoCalib and image array).",
        multiple=True,
        dimensions=(),
        deferLoad=True,
    )

    background = Input(
        name="{backgroundName}",
        storageClass="Background",
        doc="Background model for calexp.",
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

    sky_brightness_precision_table = Output(
        name="sky_brightness_precision_table",
        storageClass="ArrowAstropy",
        doc="A table containing two columns: the detector or patch IDs and the values of limiting surface "
        "brightness derived for those detectors or patches.",
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        # Update output table name for configurable dimensions
        dimen = "_visit" if "visit" in config.inputTableDimensions else "_tract"

        self.objTable = Input(
            doc=self.objTable.doc,
            name=self.objTable.name,
            storageClass=self.objTable.storageClass,
            deferLoad=self.objTable.deferLoad,
            dimensions=frozenset(sorted(config.inputTableDimensions)),
            multiple=self.objTable.multiple,
        )
        self.calexp = Input(
            doc=self.calexp.doc,
            name=self.calexp.name,
            storageClass=self.calexp.storageClass,
            deferLoad=self.calexp.deferLoad,
            dimensions=frozenset(sorted(config.inputCalibDimensions)),
            multiple=self.calexp.multiple,
        )
        self.background = Input(
            doc=self.background.doc,
            name=self.background.name,
            storageClass=self.background.storageClass,
            deferLoad=self.background.deferLoad,
            dimensions=frozenset(sorted(config.inputCalibDimensions)),
            multiple=self.background.multiple,
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
        self.sky_brightness_precision_table = Output(
            doc=self.sky_brightness_precision_table.doc,
            name=self.sky_brightness_precision_table.name + dimen,
            storageClass=self.sky_brightness_precision_table.storageClass,
            dimensions=frozenset(sorted(config.outputDataDimensions)),
        )

        assert config is not None, "Missing required config object."

        if "tract" not in config.inputTableDimensions:
            del self.skymap

        self.dimensions.update(frozenset(sorted(config.outputDataDimensions)))


class SkyBrightnessPrecisionConfig(
    PipelineTaskConfig,
    pipelineConnections=SkyBrightnessPrecisionConnections,
):
    """CConfig class for SkyBrightnessPrecisionTask"""

    inputTableDimensions = ListField[str](
        doc="Dimensions of the input object table data.",
        default=("skymap", "tract", "patch", "band"),
        optional=False,
    )
    inputCalibDimensions = ListField[str](
        doc="Dimensions of the input calibration data.",
        default=("skymap", "tract", "patch", "band"),
        optional=False,
    )
    outputDataDimensions = ListField[str](
        doc="Dimensions of the output table data.",
        default=("skymap", "tract", "patch", "band"),
        optional=False,
    )
    apertureSize = Field[int](
        doc="The size of the sky objects photometry aperture.",
        default=9,
    )
    tolerance = Field[float](
        doc="Fractional tolerance for the precision check (1% = 0.01).",
        default=0.01,
    )


class SkyBrightnessPrecisionTask(PipelineTask):
    """A task for measuring the error in the precision of the sky brightness specified by OSS-REQ-0387-V-05"""

    ConfigClass = SkyBrightnessPrecisionConfig
    _DefaultName = "skyBrightnessPrecision"

    def __init__(self, initInputs=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def _mean_in_aperture(image, x, y, aper, npix_unmasked):
        # Ensure x and y are integers and  is a numpy array
        x, y = int(x), int(y)
        img = np.asarray(image)
        rows, cols = img.shape

        # Create grid arrays for the entire image using ogrid (memory efficient)
        X, Y = np.ogrid[:rows, :cols]

        # Create the mask for points within the circle (distance from (x, y) <= aper)
        mask = (X - x) ** 2 + (Y - y) ** 2 <= aper**2

        # Compute the mean using the actual number of pixels within the mask
        nPix = np.sum(mask)

        # 9 pixel aperture has 253 pixels if unmasked, so we scale the area accordingly
        area = nPix / npix_unmasked * np.pi * aper**2

        # There are some sky sources that are out of bound, return 0 in this case
        if nPix == 0:
            meanFluxInApe = 0
        else:
            meanFluxInApe = img[mask].sum() / area

        return meanFluxInApe

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        sky_brightness_precision_struct = self.run(**{k: v for k, v in inputs.items() if k != "skymap"})
        butlerQC.put(sky_brightness_precision_struct, outputRefs)

    def run(self, objTable, calexp, background, photoCalib, wcs):
        """Output a table of per-detector or per-patch sky brightness precision measurement
        for the input visit or tract coadd

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
            Table containing per-detector or per-patch sky brightness precision values
        """

        # lookup_photoCalib = {x.dataId: x for x in photoCalib}
        # lookup_wcs = {x.dataId: x for x in wcs}
        lookup_calexp = {x.dataId: x for x in calexp}
        lookup_bg = {x.dataId: x for x in background}

        source_catalog = objTable[0].get()

        if "band" in self.config.inputCalibDimensions and len(photoCalib) > 0:
            band = photoCalib[0].dataId["band"] + "_"
        else:
            band = ""

        out = Table(
            names=["imageID", "fracWithin1pct", "maxAbsErr", "medianRatio", "nSky"],
            dtype=[np.int64, np.float64, np.float64, np.float64, np.int64],
        )

        aper = int(self.config.apertureSize)
        nAperPixels = int(SpanSet.fromShape(aper).asArray().sum())
        nPixIdeal = np.pi * aper**2

        if "visit" in self.config.inputTableDimensions:
            skyKeyName = "sky_source"
            ap_col = f"ap{aper:02d}Flux"
        else:
            skyKeyName = "sky_object"
            ap_col = f"{band}ap{aper:02d}Flux"
            
        ids_all: list[int] = []
        ratios_all: list[float] = []

        # Iterate over images (detectors or patches)
        for dataId in lookup_calexp.keys():
            idKey = "detector" if "detector" in dataId else "patch"

            # Slice catalogue rows for this image and sky sources
            isImage = source_catalog[idKey] == dataId[idKey]
            assert np.any(isImage), f"No sources found for {dataId}."
            if skyKeyName not in source_catalog.colnames:
                raise ValueError(f"No sky flag column {skyKeyName} found for {dataId}.")

            isSky = source_catalog[skyKeyName] > 0
            rows = source_catalog[isImage & isSky]

            if ap_col not in rows.colnames:
                raise ValueError(f"No aperture flux column {ap_col} found for {dataId}.")
            if len(rows) == 0:
                raise ValueError(f"No valid rows found for {dataId}.")

            cal = lookup_calexp[dataId].get()
            bg = lookup_bg[dataId].get()

            nano = cal.getPhotoCalib().instFluxToNanojansky(1.0)
            calimg = cal.image.array.astype("f8", copy=False) * nano
            x0, y0 = cal.getXY0()

            bgimg = bg.getImage().array.astype("f8", copy=False) * nano

            x = np.asarray(rows["x"], dtype="f8")
            y = np.asarray(rows["y"], dtype="f8")
            mean_flux_sky = np.asarray(rows[ap_col], dtype="f8") / nPixIdeal

            mean_img = np.empty_like(mean_flux_sky)
            mean_bg = np.empty_like(mean_flux_sky)
            for i in range(mean_img.size):
                mean_img[i] = self._mean_in_aperture(calimg, x[i] - x0, y[i] - y0, aper, nAperPixels)
                mean_bg[i] = self._mean_in_aperture(bgimg, x[i] - x0, y[i] - y0, aper, nAperPixels)

            good = mean_bg != 0.0
            if not np.any(good):
                raise ValueError(f"No valid background flux found for {dataId}.")

            # SBRatio: (background + skyFlux) / background
            sb_ratio = (mean_bg[good] + mean_flux_sky[good]) / mean_bg[good]
            img_id = int(dataId[idKey])
            ids_all.extend([img_id] * sb_ratio.size)
            ratios_all.extend(sb_ratio.astype("f8").tolist())

        out["imageID"] = np.array(ids_all, dtype=np.int64)
        out["sb_ratio"] = np.array(ratios_all, dtype=np.float64)
            

            # tol = float(self.config.tolerance)
            # fracWithin = np.count_nonzero((sb_ratio >= (1 - tol)) & (sb_ratio <= (1 + tol))) / sb_ratio.size
            # maxAbsErr = np.max(np.abs(sb_ratio - 1.0))
            # medianRatio = np.median(sb_ratio)

            # out.add_row(
            #     [
            #         int(dataId[idKey]),
            #         float(fracWithin),
            #         float(maxAbsErr),
            #         float(medianRatio),
            #         int(sb_ratio.size),
            #     ]
            # )

        return Struct(sky_brightness_precision_table=out)


class SkyBrightnessPrecisionAnalysisConnections(
    AnalysisBaseConnections,
    defaultTemplates={"outputName": "skyBrightnessPrecision"},
):
    data = Input(
        name="sky_brightness_precision_table",
        storageClass="ArrowAstropy",
        doc="A table containing sky brightness precision metrics.",
        deferLoad=True,
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        dimen = "_visit" if "visit" in config.inputDataDimensions else "_tract"

        self.data = Input(
            name=self.data.name + dimen,
            storageClass=self.data.storageClass,
            doc=self.data.doc,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputDataDimensions)),
        )

        self.dimensions.update(frozenset(sorted(config.inputDataDimensions)))


class SkyBrightnessPrecisionAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=SkyBrightnessPrecisionAnalysisConnections
):
    inputDataDimensions = ListField[str](
        doc="Dimensions of the input table data.",
        default=(),
        optional=False,
    )


class SkyBrightnessPrecisionAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SkyBrightnessPrecisionAnalysisConfig
    _DefaultName = "skyBrightnessPrecisionAnalysis"
