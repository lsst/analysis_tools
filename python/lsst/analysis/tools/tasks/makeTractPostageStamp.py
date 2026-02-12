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
    "MakeTractPostageStampConfig",
    "MakeTractPostageStampTask",
)

import copy

import matplotlib.cm as cm
import matplotlib.patheffects as pathEffects
import numpy as np

import lsst.pipe.base as pipeBase
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap
from lsst.utils.plotting import make_figure

from ..utils import getPatchCorners, getTractCorners


class MakeTractPostageStampConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("skymap", "tract", "band"),
    defaultTemplates={
        "inputName": "deepCoadd_calexpBin",
        "outputName": "deepTract_PostageStamp",
    },
):
    data = ct.Input(
        doc="Binned deepCoadd calibrated exposures to read from the butler.",
        name="{inputName}",
        storageClass="ExposureF",
        deferLoad=True,
        dimensions=(
            "skymap",
            "tract",
            "patch",
            "band",
        ),
        multiple=True,
    )

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    postageStamp = ct.Output(
        doc="A postagestamp composite image of the tract.",
        name="{outputName}",
        storageClass="Plot",
        dimensions=("skymap", "tract"),
    )


class MakeTractPostageStampConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=MakeTractPostageStampConnections
):
    pass


class MakeTractPostageStampTask(pipeBase.PipelineTask):

    ConfigClass = MakeTractPostageStampConfig
    _DefaultName = "makeTractPostageStamp"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Takes a set of coadded patch Exposures and displays them
        in their corresponding positions within a tract.
        Empty patches - those that do not have any coverage - are shown
        as hatched squares.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
            A butler which is specialized to operate in the context of a
            `lsst.daf.butler.Quantum`.
        inputRefs : `lsst.pipe.base.InputQuantizedConnection`
            Datastructure containing named attributes 'data and 'skymap'.
            The values of these attributes are the corresponding
            `lsst.daf.butler.DatasetRef` objects defined in the corresponding
            `PipelineTaskConnections` class.
        outputRefs : `lsst.pipe.base.OutputQuantizedConnection`
            Datastructure containing named attribute 'postageStamp'.
            The value of this attribute is the corresponding
            `lsst.daf.butler.DatasetRef` object defined in the corresponding
            `PipelineTaskConnections` class.
        """

        inputs = butlerQC.get(inputRefs)
        patches = inputs["data"]
        skymap = inputs["skymap"]
        band = butlerQC.quantum.dataId["band"]
        tract = butlerQC.quantum.dataId["tract"]

        fig = self.makeTractPostageStamp(skymap, tract, patches, band)

        butlerQC.put(pipeBase.Struct(postageStamp=fig), outputRefs)

    def makeTractPostageStamp(self, skymap, tract, patches, band):
        """Takes the coadd patch exposures and displays them on a
        set of axes. The axes boundaries are those of the tract.
        Patches are annoted with their patch identification number.
        Empty patches - those that do not have any coverage - are shown
        as hatched squares.

        Parameters
        ----------
        skymap : `lsst.skymap`
        tract : `int`
            Identification number of tract.
        patches : `list` [`DeferredDatasetHandle`]
            List of handles for patch coadd exposures to display.
        band : `str`
            Filter band. Only used to annotate the plot.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            Plot displaying the tract with coadd exposures displayed.
        """

        tractInfo = skymap.generateTract(tract)
        tractCorners = getTractCorners(skymap, tract)
        tractRas = [ra for (ra, dec) in tractCorners]
        RaSpansZero = max(tractRas) > 360.0

        cmap = cm.grey
        cmap.set_bad(alpha=0)
        cmapred = cm.Reds
        cmapred.set_bad(alpha=0)

        fig = make_figure(dpi=300)
        ax = fig.add_subplot(111)

        emptyPatches = np.arange(tractInfo.getNumPatches()[0] * tractInfo.getNumPatches()[1]).tolist()

        for patch in patches:

            patchId = patch.dataId["patch"]
            emptyPatches.remove(patchId)

            # Create the patch axes at the appropriate location in tract:
            patchCorners = getPatchCorners(tractInfo, patchId)
            ras = [ra for (ra, dec) in patchCorners]
            decs = [dec for (ra, dec) in patchCorners]

            # Account for the RA wrapping using negative RA values.
            # This is rectified when the final axes are built.
            if RaSpansZero:
                ras = [ra - 360 if ra > 180.0 else ra for ra in ras]
            Extent = (max(ras), min(ras), max(decs), min(decs))
            ax.plot(
                [min(ras), max(ras), max(ras), min(ras), min(ras)],
                [min(decs), min(decs), max(decs), max(decs), min(decs)],
                "k",
                lw=0.5,
            )

            # Fetch the images and display:
            im = patch.get(component="image").array
            mask = patch.get(component="mask").array & 2**8 > 0

            imdata = copy.deepcopy(im)
            imdata[mask] = np.nan
            med = np.nanmedian(imdata)
            mad = np.nanmedian(np.fabs(im - med))
            vmin = med - 1 * mad
            vmax = med + 15 * mad
            ax.imshow(imdata, extent=Extent, vmin=vmin, vmax=vmax, cmap=cmap)

            # Highlight regions of "NODATA".
            imnodata = copy.deepcopy(im)
            imnodata[~mask] = np.nan
            ax.imshow(imnodata, extent=Extent, alpha=0.8, cmap=cmapred)

            ax.annotate(
                patchId,
                (np.mean(ras), np.mean(decs)),
                color="k",
                ha="center",
                va="center",
                path_effects=[pathEffects.withStroke(linewidth=2, foreground="w")],
            )

        # Indicate blank patches as hatched regions
        for patchId in emptyPatches:

            patchCorners = getPatchCorners(tractInfo, patchId)
            ras = [ra for (ra, dec) in patchCorners]
            decs = [dec for (ra, dec) in patchCorners]

            # Account for the RA wrapping using negative RA values.
            if RaSpansZero:
                ras = [ra - 360 if ra > 180.0 else ra for ra in ras]

            Extent = (max(ras), min(ras), max(decs), min(decs))
            ax.plot(
                [min(ras), max(ras), max(ras), min(ras), min(ras)],
                [min(decs), min(decs), max(decs), max(decs), min(decs)],
                "k",
                lw=0.5,
            )

            cs = ax.contourf(np.ones((10, 10)), 1, hatches=["xx"], extent=Extent, colors="none")
            cs.set_edgecolors("red")
            ax.annotate(
                patchId,
                (np.mean(ras), np.mean(decs)),
                color="k",
                ha="center",
                va="center",
                path_effects=[pathEffects.withStroke(linewidth=2, foreground="w")],
            )

        # Draw axes around the entire tract:
        ax.set_xlabel("R.A. (deg)", fontsize=20)
        ax.set_ylabel("Dec. (deg)", fontsize=20)

        tractRas = [ra for (ra, dec) in tractCorners]
        # Account for the RA wrapping using negative RA values.
        if RaSpansZero:
            tractRas = [ra - 360.0 for ra in tractRas]

        ax.set_xlim(max(tractRas), min(tractRas))
        ticks = [t for t in ax.get_xticks() if t >= min(tractRas) and t <= max(tractRas)]

        # Rectify potential negative RA values via tick labels
        tickLabels = [f"{t % 360:.1f}" for t in ticks]
        ax.set_xticks(ticks, tickLabels)

        tractDecs = [dec for (ra, dec) in tractCorners]
        ax.set_ylim(min(tractDecs), max(tractDecs))

        ax.tick_params(axis="both", labelsize=10, length=0, pad=1.5)
        ax.set_title(str(tract) + ": " + band, fontsize=20)
        fig.canvas.draw()

        return fig
