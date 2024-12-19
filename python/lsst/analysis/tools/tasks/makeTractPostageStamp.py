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

import lsst.pipe.base as pipeBase
import matplotlib.cm as cm
import matplotlib.patheffects as pathEffects
import numpy as np
from lsst.geom import Box2D
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap
from lsst.utils.plotting import make_figure


class MakeTractPostageStampConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("skymap", "tract", "band"),
    defaultTemplates={
        "inputName": "deepCoadd_calexpBin",  # Do we want to restrict to "deep" coadds?
        "outputName": "deepTract_PostageStamp",
    },  # Do we want "Tract" in the name?
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

        inputs = butlerQC.get(inputRefs)
        patches = inputs["data"]
        skymap = inputs["skymap"]
        band = butlerQC.quantum.dataId["band"]
        tract = butlerQC.quantum.dataId["tract"]

        tractInfo = skymap.generateTract(tract)

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
            pRas, pDecs = self.getPatchCorners(tractInfo, patchId)
            pExtent = (max(pRas), min(pRas), max(pDecs), min(pDecs))
            ax.plot(pRas + [pRas[0]], pDecs + [pDecs[0]], "k", lw=0.5)

            # Fetch the images and display:
            im = patch.get(component="image").array
            mask = patch.get(component="mask").array & 2**8 > 0

            imdata = copy.deepcopy(im)
            imdata[mask] = np.nan
            med = np.nanmedian(imdata)
            mad = np.nanmedian(np.fabs(im - med))
            vmin = med - 1 * mad
            vmax = med + 15 * mad
            ax.imshow(imdata, extent=pExtent, vmin=vmin, vmax=vmax, cmap=cmap)

            # Highlight regions of "NODATA".
            imnodata = copy.deepcopy(im)
            imnodata[~mask] = np.nan
            ax.imshow(imnodata, extent=pExtent, alpha=0.8, cmap=cmapred)

            ax.annotate(
                patchId,
                (np.mean(pRas), np.mean(pDecs)),
                color="k",
                ha="center",
                va="center",
                path_effects=[pathEffects.withStroke(linewidth=2, foreground="w")],
            )

        # Indicate blank patches as hatched regions
        for patchId in emptyPatches:

            pRas, pDecs = self.getPatchCorners(tractInfo, patchId)
            pExtent = (max(pRas), min(pRas), max(pDecs), min(pDecs))
            ax.plot(pRas + [pRas[0]], pDecs + [pDecs[0]], "k", lw=0.5)

            cs = ax.contourf(np.ones((10, 10)), 1, hatches=["xx"], extent=pExtent, colors="none")
            for c in cs.collections:
                c.set_edgecolor("red")
            ax.annotate(
                patchId,
                (np.mean(pRas), np.mean(pDecs)),
                color="k",
                ha="center",
                va="center",
                path_effects=[pathEffects.withStroke(linewidth=2, foreground="w")],
            )

        # Draw axes around the entire tract:
        tCorners = self.getTractCorners(skymap, tractInfo.getId())
        tRas = [ra for (ra, dec) in tCorners]
        tDecs = [dec for (ra, dec) in tCorners]

        ax.set_xlabel("R.A. (deg)", fontsize=20)
        ax.set_ylabel("Dec. (deg)", fontsize=20)
        ax.tick_params(axis="both", labelsize=10, length=0, pad=1.5)
        ax.set_xlim(max(tRas), min(tRas))
        ax.set_ylim(min(tDecs), max(tDecs))
        ax.set_title(str(tract) + ": " + band, fontsize=20)
        fig.canvas.draw()

        butlerQC.put(pipeBase.Struct(postageStamp=fig), outputRefs)

    def getTractCorners(self, skymap, tract):
        """Calculate the corners of a tract, given  skymap.
        ***This needs to be moved from here***

        Parameters
        ----------
        skymap : `lsst.skymap`
        tract : `int`

        Returns
        -------
        corners : `list` of `tuples` of `float`

        Notes
        -----
        Corners are returned in degrees and wrapped in ra.
        """
        # Find the tract corners
        tractCorners = skymap[tract].getVertexList()
        corners = [(corner.getRa().asDegrees(), corner.getDec().asDegrees()) for corner in tractCorners]
        minRa = np.min([corner[0] for corner in corners])
        maxRa = np.max([corner[0] for corner in corners])

        # If the tract needs wrapping in ra, wrap it
        if maxRa - minRa > 10:
            x = maxRa
            maxRa = 360 + minRa
            minRa = x
            minDec = np.min([corner[1] for corner in corners])
            maxDec = np.max([corner[1] for corner in corners])
            corners = [(minRa, minDec), (maxRa, minDec), (maxRa, maxDec), (minRa, maxDec)]
        return corners

    def getPatchCorners(self, tractInfo, patchId):

        patchInfo = tractInfo.getPatchInfo(patchId)
        pcorners = Box2D(patchInfo.getInnerBBox()).getCorners()
        tractWcs = tractInfo.getWcs()

        skyCoords = tractWcs.pixelToSky(pcorners)

        ras = [ra.asDegrees() for (ra, dec) in skyCoords]
        decs = [dec.asDegrees() for (ra, dec) in skyCoords]

        return (ras, decs)
