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
    "CellCoaddPsfSurveyConfig",
    "CellCoaddPsfSurveyTask",
)

import numpy as np
from astropy.table import Table

import lsst.geom as geom
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base import connectionTypes as cT
from lsst.pipe.tasks.coaddBase import makeSkyInfo
from lsst.skymap import BaseSkyMap


class CellCoaddPsfSurveyConnections(
    PipelineTaskConnections,
    dimensions=("tract", "band", "skymap"),
    defaultTemplates={"coaddName": "deep", "coaddSuffix": "Cell"},
):
    coadd = cT.Input(
        doc="Cell coadd exposure from which PSF and mask are read.",
        name="{coaddName}_coadd",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "band", "skymap"),
        multiple=True,
        deferLoad=True,
    )
    psfWeightFraction = cT.Input(
        doc="Per-cell PSF weight fraction image (temporary diagnostic output from "
            "AssembleCellCoaddTask). Absent when do_psf_weight_fraction_image is False.",
        name="{coaddName}Coadd{coaddSuffix}_psfWeightFraction",
        storageClass="ImageF",
        dimensions=("tract", "patch", "band", "skymap"),
        multiple=True,
        deferLoad=True,
        minimum=0,
    )
    skyMap = cT.PrerequisiteInput(
        doc="Sky map used to locate cell centres within each patch.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )
    psfCellMetrics = cT.Output(
        doc="Per-cell table of PSF metrics for the whole tract.",
        name="{coaddName}_coadd_cell_psf_metrics",
        storageClass="ArrowAstropy",
        dimensions=("tract", "band", "skymap"),
    )


class CellCoaddPsfSurveyConfig(PipelineTaskConfig, pipelineConnections=CellCoaddPsfSurveyConnections):
    pass


class CellCoaddPsfSurveyTask(PipelineTask):
    """Collect per-cell PSF metrics from a tract of cell coadds.

    For every cell in every patch of the tract the task evaluates the PSF
    determinant radius (psfSigma) at the inner-bbox centre, records the
    inner-bbox geometry, reads the INEXACT_PSF mask flag at that centre pixel,
    and (when available) reads the PSF weight fraction produced by
    AssembleCellCoaddTask.  All cells from all patches are consolidated into a
    single tract-level ArrowAstropy table.

    Parameters
    ----------
    Columns in the output table:

    ``patch``
        Sequential patch index (int).
    ``cell_ix``, ``cell_iy``
        Cell grid indices within the patch (int).  Together with ``patch``
        these uniquely identify a cell and are sufficient to recover its
        geometry from the skymap.
    ``psf_ixx``, ``psf_iyy``, ``psf_ixy``
        Second moments of the PSF evaluated at the cell centre (float, pixels²).
        NaN when PSF evaluation fails.
    ``psf_sigma``
        PSF determinant radius ``(Ixx*Iyy - Ixy²)^0.25`` in pixels (float).
        NaN when PSF evaluation fails.
    ``pwf``
        Mean PSF weight fraction over the inner bbox (float).
        NaN when the psfWeightFraction dataset is absent.
    ``inexact_psf``
        Whether the INEXACT_PSF mask bit is set at the centre pixel (bool).
    """

    ConfigClass = CellCoaddPsfSurveyConfig
    _DefaultName = "cellCoaddPsfSurvey"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        sky_map = butlerQC.get(inputRefs.skyMap)
        coadd_handles = butlerQC.get(inputRefs.coadd)
        pwf_handles = butlerQC.get(inputRefs.psfWeightFraction)

        pwf_by_patch = {int(h.dataId["patch"]): h for h in pwf_handles}

        tract_id = butlerQC.quantum.dataId["tract"]
        result = self.run(coadd_handles, pwf_by_patch, sky_map, tract_id)
        butlerQC.put(result, outputRefs)

    def run(self, coadd_handles, pwf_by_patch, sky_map, tract_id):
        """Build the per-cell PSF metrics table.

        Parameters
        ----------
        coadd_handles : `list` of `DeferredDatasetHandle`
            Deferred handles for the per-patch ``deep_coadd`` ExposureF datasets.
        pwf_by_patch : `dict` [`int`, `DeferredDatasetHandle`]
            Deferred handles for the per-patch ``psfWeightFraction`` ImageF
            datasets, keyed by sequential patch index.  Empty when the dataset
            is not available.
        sky_map : `BaseSkyMap`
            Sky map providing patch and cell geometry.
        tract_id : `int`
            Tract identifier used to look up patch info in the sky map.

        Returns
        -------
        result : `Struct`
            ``psfCellMetrics`` — an `astropy.table.Table` with one row per cell.
        """
        rows = []
        for coadd_handle in coadd_handles:
            patch_id = int(coadd_handle.dataId["patch"])
            exp = coadd_handle.get()
            psf = exp.getPsf()
            mask = exp.getMask()

            try:
                inexact_psf_bit = mask.getPlaneBitMask("INEXACT_PSF")
            except Exception:
                inexact_psf_bit = 0

            pwf_image = pwf_by_patch[patch_id].get() if patch_id in pwf_by_patch else None

            sky_info = makeSkyInfo(sky_map, tractId=tract_id, patchId=patch_id)

            for cell_info in sky_info.patchInfo:
                inner = cell_info.inner_bbox
                center = geom.Point2D(inner.getCenter())
                center_i = geom.Point2I(center)

                try:
                    shape = psf.computeShape(center)
                    psf_ixx = float(shape.getIxx())
                    psf_iyy = float(shape.getIyy())
                    psf_ixy = float(shape.getIxy())
                    psf_sigma = float(shape.getDeterminantRadius())
                except Exception:
                    psf_ixx = psf_iyy = psf_ixy = psf_sigma = np.nan

                pwf = np.nan
                if pwf_image is not None:
                    inner_arr = pwf_image[inner].array
                    finite = inner_arr[np.isfinite(inner_arr)]
                    if finite.size:
                        pwf = float(finite.mean())

                inexact_psf = bool(int(mask[center_i]) & inexact_psf_bit)

                rows.append(
                    {
                        "patch": patch_id,
                        "cell_ix": cell_info.index.x,
                        "cell_iy": cell_info.index.y,
                        "psf_ixx": psf_ixx,
                        "psf_iyy": psf_iyy,
                        "psf_ixy": psf_ixy,
                        "psf_sigma": psf_sigma,
                        "pwf": pwf,
                        "inexact_psf": inexact_psf,
                    }
                )

        if not rows:
            table = Table(
                names=["patch", "cell_ix", "cell_iy",
                       "psf_ixx", "psf_iyy", "psf_ixy", "psf_sigma", "pwf", "inexact_psf"],
                dtype=["i4", "i4", "i4", "f8", "f8", "f8", "f8", "f8", "bool"],
            )
        else:
            table = Table(rows=rows)

        return Struct(psfCellMetrics=table)
