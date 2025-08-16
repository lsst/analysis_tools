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
    dimensions=(),
    defaultTamplates={
        "detectionTableName": "object_all",
        "calexpName": "deep_coadd.calexp",
        "backgroundName": "deep_coadd.calexpBackground",
        "photoCalibName": "deep_coadd.photoCalib",
        "wcsName": "deep_coadd.wcs",
    },
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

    calexpBackground = Input(
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
        "brightness derived for those detectors or patches.",  # TODO: change doc
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        # Update output table name for configurable dimensions
        dimen = "_visit" if "visit" in config.inputTableDimensions else "_tract"
        self.objTable = Input(
            doc=self.data.doc,
            name=self.data.name,
            storageClass=self.data.storageClass,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputTableDimensions)),
            multiple=self.data.multiple,
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
    def _disk_mask(radius_px: int):
        r = int(radius_px)
        yy, xx = np.ogrid[-r : r + 1, -r : r + 1]
        return (xx * xx + yy * yy) <= r * r

    @staticmethod
    def _mean_in_aperture(
        image_array: np.ndarray, x: float, y: float, disk_mask: np.ndarray, nAperPixels: int, aper: int
    ) -> float:
        """Mean over an *effective* area for circular aperture centered at (x,y)."""
        r = disk_mask.shape[0] // 2
        x0 = int(np.round(x)) - r
        y0 = int(np.round(y)) - r
        x1 = x0 + disk_mask.shape[1]
        y1 = y0 + disk_mask.shape[0]

        H, W = image_array.shape
        xs0, ys0 = max(0, x0), max(0, y0)
        xs1, ys1 = min(W, x1), min(H, y1)
        if xs0 >= xs1 or ys0 >= ys1:
            return 0.0

        im_cut = image_array[ys0:ys1, xs0:xs1]
        mk_cut = disk_mask[
            (ys0 - y0) : (ys1 - y1 + disk_mask.shape[0]), (xs0 - x0) : (xs1 - x1 + disk_mask.shape[1])
        ]
        mk_cut = mk_cut[: im_cut.shape[0], : im_cut.shape[1]]

        nPix = int(mk_cut.sum())
        if nPix == 0:
            return 0.0

        ideal_area = np.pi * (aper**2)
        eff_area = ideal_area * (nPix / nAperPixels)
        return float(im_cut[mk_cut].sum() / eff_area)

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

            if "visit" in self.config.inputTableDimensions:
                skyKeyName = "sky_source"
            else:
                skyKeyName = "sky_object"

            out = Table(
                names=["imageID", "fracWithin1pct", "maxAbsErr", "medianRatio", "nSky"],
                dtype=[np.int64, np.float64, np.float64, np.float64, np.int64],
            )

            aper = int(self.config.apertureSize)
            ap_col = f"{band}ap{aper:02d}Flux"
            disk = self._disk_mask(aper)
            nAperPixels = int(SpanSet.fromShape(aper).asArray().sum())
            nPixIdeal = np.pi * aper**2

            # Iterate over images (detectors or patches)
            for dataId in lookup_calexp.keys():
                idKey = "detector" if "detector" in dataId else "patch"

                # Slice catalogue rows for this image and sky sources
                isImage = source_catalog[idKey] == dataId[idKey]
                if skyKeyName not in source_catalog.colnames:
                    raise ValueError(f"No sky flag column {skyKeyName} found for {dataId}.")

                isSky = source_catalog[skyKey] > 0
                rows = source_catalog[isImage & isSky]

                if ap_col not in rows.colnames:
                    raise ValueError(f"No aperture flux column {ap_col} found for {dataId}.")
                if len(rows) == 0:
                    out.add_row([int(dataId[idKey]), np.nan, np.nan, np.nan, 0])
                    continue

                cal = lookup_calexp[dataId].get()
                bg = lookup_bg[dataId].get()

                nano = cal.getPhotoCalib().instFluxToNanojansky(1.0)
                calimg = cal.image.array.astype("f8", copy=False) * nano
                bgimg = bg.getImage().array.astype("f8", copy=False) * nano

                x = np.asarray(rows["x"], dtype="f8")
                y = np.asarray(rows["y"], dtype="f8")
                mean_flux_sky = np.asarray(rows[ap_col], dtype="f8") / nPixIdeal

                mean_img = np.empty_like(mean_flux_sky)
                mean_bg = np.empty_like(mean_flux_sky)
                for i in range(mean_img.size):
                    mean_img[i] = self._mean_in_aperture(calimg, x[i], y[i], disk, nAperPixels, aper)
                    mean_bg[i] = self._mean_in_aperture(bgimg, x[i], y[i], disk, nAperPixels, aper)

                good = mean_bg != 0.0
                if not np.any(good):
                    out.add_row([int(dataId[idKey]), np.nan, np.nan, np.nan, int(mean_bg.size)])
                    continue

                # SBRatio per the notebook: (background + skyFlux) / background
                sb_ratio = (mean_bg[good] + mean_flux_sky[good]) / mean_bg[good]

                tol = float(self.config.tolerance)
                fracWithin = (
                    np.count_nonzero((sb_ratio >= (1 - tol)) & (sb_ratio <= (1 + tol))) / sb_ratio.size
                )
                maxAbsErr = np.max(np.abs(sb_ratio - 1.0))
                medianRatio = np.median(sb_ratio)

                out.add_row(
                    [
                        int(dataId[idKey]),
                        float(fracWithin),
                        float(maxAbsErr),
                        float(medianRatio),
                        int(sb_ratio.size),
                    ]
                )

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
