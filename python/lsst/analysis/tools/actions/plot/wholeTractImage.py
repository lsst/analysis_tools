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

__all__ = ("WholeTractImage",)

import logging
from typing import Mapping, Optional

import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ImageNormalize
from lsst.pex.config import (
    ChoiceField,
    Field,
    FieldValidationError,
    ListField,
)
from lsst.pex.config.configurableActions import ConfigurableActionField
from lsst.skymap import BaseSkyMap
from lsst.utils.plotting import make_figure, set_rubin_plotstyle
from matplotlib.figure import Figure

from ...interfaces import (
    KeyedData,
    KeyedDataSchema,
    PlotAction,
    TensorAction,
    VectorAction,
)
from ...utils import getPatchCorners, getTractCorners
from .calculateRange import Asinh, Perc

_LOG = logging.getLogger(__name__)


class WholeTractImage(PlotAction):
    """
    Produces a figure displaying whole-tract coadd pixel data as a 2D image.

    The figure is constructed from all patches covering the tract. Regions of
    NO_DATA or where no coadd exists are shown as red shading or red hatches,
    respectively.

    Either the image, pixel mask, or variance components of the coadd can be
    displayed. In the case of the pixel mask, one or more bitmaskPlanes must
    be specified; the specified bitmaskPlanes are OR-combined, with flagged
    pixels given a value of 1, and unflagged pixels given a value of 1.
    """

    component = ChoiceField[str](
        doc="Coadd component to display. Can take one of image, mask, variance. Default: image.",
        default="image",
        allowed={plane: plane for plane in ("image", "mask", "variance")},
    )

    bitmaskPlanes = ListField[str](
        doc="List of names of bitmask plane(s) to display when displaying the "
        "mask plane. Bitmask planes are OR-combined. Flagged pixels are given "
        "a value of 1; unflagged pixels are given a value of 0. "
        "Optional when displaying either the image or variance planes. "
        "Required when displaying the mask plane.",
        optional=True,
    )

    showPatchIds = Field[bool](
        doc="Show the patch IDs in the centre of each patch. Default: False",
        default=False,
    )

    showColorbar = Field[bool](
        doc="Show a colorbar alongside the main plot. Default: False",
        default=False,
    )

    zAxisLabel = Field[str](
        doc="Label to display on the colorbar. Optional",
        optional=True,
    )

    interval = ConfigurableActionField[VectorAction](
        doc="Action to calculate the min and max values of the image scale. Default: Perc.",
        default=Perc,
    )

    colorbarCmap = ChoiceField[str](
        doc="Matplotlib colormap to use for the displayed image. Default: gray",
        default="gray",
        allowed={name: name for name in plt.colormaps()},
    )

    noDataColor = Field[str](
        doc="Matplotlib color to use to indicate regions of no data. Default: red",
        default="red",
    )

    noDataValue = Field[int](
        doc="If data doesn't contain a mask plane, the value in the image plane to "
        "assign the noDataColor to. Optional.",
        optional=True,
    )

    vmaxFloor = Field[float](
        doc="The floor of the vmax value of the colorbar",
        default=None,
        optional=True,
    )

    stretch = ConfigurableActionField[TensorAction](
        doc="Action to calculate the stretch of the image scale. Default: Asinh",
        default=Asinh,
    )

    displayAsPostageStamp = Field[bool](
        doc="Display as a figure to be used as postage stamp. No plotInfo or legend is shown, "
        "and large fonts are used for axis labels.",
        default=False,
    )

    def validate(self):
        super().validate()

        if self.component == "mask" and self.bitmaskPlanes is None:
            raise FieldValidationError(
                self.__class__.bitmaskPlanes,
                self,
                "'bitmaskPlanes' must be specified if displaying the mask plane.",
            )
        if self.bitmaskPlanes is not None and self.component != "mask":
            raise FieldValidationError(
                self.__class__.component,
                self,
                "'component' must be set to the mask plane if 'bitmaskPlanes' is specified.",
            )

    def getInputSchema(self) -> KeyedDataSchema:
        base = []
        base.append((self.component, KeyedData))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makeFigure(data, **kwargs)

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        needed = self.getInputSchema()
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in the input data")

    def makeFigure(
        self,
        data: KeyedData,
        tractId: int,
        skymap: BaseSkyMap,
        plotInfo: Optional[Mapping[str, str]] = None,
        **kwargs,
    ) -> Figure:
        """Make a figure displaying the input pixel data.

        Parameters
        ----------
        data : `lsst.analysis.tools.interfaces.KeyedData`
            A python dict-of-dicts containing the pixel data to display in the
            figure. The top level keys are named after the coadd component(s),
            and must contain at least 'mask'. The next level keys are named
            after the patch ID of the coadd component contained as their
            corresponding value.
        tractId : `int`
            Identification number of the tract to be displayed.
        skymap : `lsst.skymap.BaseSkyMap`
            The sky map used for this dataset. This is referred-to to determine
            the location of the tract on-sky (for RA and Dec axis ranges) and
            the location of the patches within the tract.
        plotInfo : `dict`, optional
            A dictionary of information about the data being plotted with keys:

            ``"run"``
                The output run for the plots (`str`).
            ``"skymap"``
                The type of skymap used for the data (`str`).
            ``"band"``
                The filter used for this data (`str`). Optional
            ``"tract"``
                The tract that the data comes from (`str`).

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.

        Examples
        --------
        An example wholeTractImage plot may be seen below:

        .. image:: /_static/analysis_tools/wholeTractImageExample.png

        For further details on how to generate a plot, please refer to the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """

        tractInfo = skymap.generateTract(tractId)
        tractCorners = getTractCorners(skymap, tractId)
        tractRas = [ra for (ra, dec) in tractCorners]
        RaSpansZero = max(tractRas) > 360.0

        cmap = cm.get_cmap(self.colorbarCmap).reversed().copy()
        cmap.set_bad(self.noDataColor, alpha=0.6 if self.noDataColor == "red" else 1.0)

        set_rubin_plotstyle()
        fig = make_figure()
        ax = fig.add_subplot(111)

        if plotInfo is None:
            plotInfo = {}
        plotInfo["component"] = self.component

        if self.bitmaskPlanes is not None:
            plotInfo["maskPlanes"] = self.bitmaskPlanes

        if self.displayAsPostageStamp:
            axisLabelFontSize = 20
            tickMarkFontSize = 10
            boundaryColor = "k"
            boundaryAlpha = 1.0
            boundaryWidth = 0.5
        else:
            axisLabelFontSize = 8
            tickMarkFontSize = 8
            boundaryColor = "r" if "viridis" in self.colorbarCmap.lower() else "c"
            boundaryAlpha = 0.3
            boundaryWidth = 1.0

        # Keep a record of the "empty" patches that do not have coadds.
        emptyPatches = np.arange(tractInfo.getNumPatches()[0] * tractInfo.getNumPatches()[1]).tolist()

        # Extract the pixel arrays for all patches prior to plotting.
        # This allows for a global image normalisation to be calculated.
        imStack = dict()
        allPix = np.array([])
        patchIds = data[self.component].keys()
        first = True
        for patchId in patchIds:

            if first:
                if "mask" in data:
                    noDataBitmask = data["mask"][patchId].getPlaneBitMask("NO_DATA")
                    maskPlanes = set(data["mask"][patchId].getMaskPlaneDict())
                    bitmaskPlanes = set(self.bitmaskPlanes) if self.bitmaskPlanes else set()
                    if bitmaskPlanes:
                        if missingMaskPlanes := bitmaskPlanes - maskPlanes:
                            _LOG.info(
                                "%s not found among the mask planes for patchId=%d",
                                missingMaskPlanes,
                                patchId,
                            )
                        else:
                            first = False
                        bitmasks = data["mask"][patchId].getPlaneBitMask(bitmaskPlanes & maskPlanes)

            emptyPatches.remove(patchId)
            im = data[self.component][patchId].array
            if self.bitmaskPlanes:
                im = (im & bitmasks > 0) * 1.0

            if "mask" in data:
                noDataMask = data["mask"][patchId].array & noDataBitmask > 0
            elif self.noDataValue is not None:
                noDataMask = data[self.component][patchId].array == self.noDataValue
            else:
                noDataMask = np.zeros_like(data[self.component][patchId].array) > 0

            allPix = np.append(allPix, im[~noDataMask].flatten())
            imStack[patchId] = np.ma.masked_array(im, mask=noDataMask)

        # It is possible that all pixels are flagged NO_DATA.
        # In which case, set vmin & vmax to arbitrary values.
        if len(allPix) == 0:
            vmin, vmax = (0, 1)
        else:
            vmin, vmax = self.interval(allPix)

        # Set a floor to vmax. Useful for low dymanic range data.
        if self.vmaxFloor is not None:
            vmax = max(vmax, self.vmaxFloor)

        for patchId, im in imStack.items():

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
                boundaryColor,
                lw=boundaryWidth,
                alpha=boundaryAlpha,
            )

            norm = ImageNormalize(vmin=vmin, vmax=vmax)
            stretchedIm = self.stretch(norm(im))
            masked_stretched = np.ma.masked_array(
                norm.inverse(stretchedIm.data),
                mask=stretchedIm.mask,
            )
            plotIm = ax.imshow(masked_stretched, vmin=vmin, vmax=vmax, extent=Extent, cmap=cmap)

            if self.showPatchIds:
                ax.annotate(
                    patchId,
                    (np.mean(ras), np.mean(decs)),
                    color="k",
                    ha="center",
                    va="center",
                    path_effects=[pathEffects.withStroke(linewidth=2, foreground="w")],
                )

        # Indicate the empty patches with red hatching
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
                boundaryColor,
                lw=boundaryWidth,
                alpha=boundaryAlpha,
            )

            cs = ax.contourf(np.ones((10, 10)), 1, hatches=["xx"], extent=Extent, colors="none")
            cs.set_edgecolors("red")
            if self.showPatchIds:
                ax.annotate(
                    patchId,
                    (np.mean(ras), np.mean(decs)),
                    color="k",
                    ha="center",
                    va="center",
                    path_effects=[pathEffects.withStroke(linewidth=2, foreground="w")],
                )

        # Draw axes around the entire tract:
        ax.set_xlabel("R.A. (deg)", fontsize=axisLabelFontSize)
        ax.set_ylabel("Dec. (deg)", fontsize=axisLabelFontSize)

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

        ax.tick_params(axis="both", labelsize=tickMarkFontSize, length=0, pad=1.5)

        if self.showColorbar:
            cax = fig.add_axes([0.90, 0.11, 0.04, 0.77])
            cbar = fig.colorbar(plotIm, cax=cax, extend="both")
            cbar.ax.tick_params(labelsize=tickMarkFontSize)
            if self.zAxisLabel:
                colorbarLabel = self.zAxisLabel
            else:
                colorbarLabel = ""
            text = cax.text(
                0.5,
                0.5,
                colorbarLabel,
                color="k",
                rotation="vertical",
                transform=cax.transAxes,
                ha="center",
                va="center",
                fontsize=10,
            )
            text.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])

        if self.displayAsPostageStamp:
            if "band" in plotInfo:
                title = f"{str(tractId)}; {plotInfo['band']}"
            else:
                title = f"{str(tractId)}"
            ax.set_title(title, fontsize=20)

        if not self.displayAsPostageStamp:
            if "mask" in data:
                noDataPatch = patches.Rectangle(
                    (0.8, 1.1), 0.05, 0.04, transform=ax.transAxes, facecolor="red", alpha=0.6, clip_on=False
                )
                ax.add_patch(noDataPatch)
                ax.text(0.86, 1.115, "NO_DATA", transform=ax.transAxes, va="center", ha="left", fontsize=8)

            noCoaddPatch = patches.Rectangle(
                (0.8, 1.02),
                0.05,
                0.04,
                transform=ax.transAxes,
                facecolor="none",
                edgecolor="red",
                hatch="xx",
                clip_on=False,
            )
            ax.add_patch(noCoaddPatch)
            ax.text(0.86, 1.04, "No coadd", transform=ax.transAxes, va="center", ha="left", fontsize=8)

            fig = addPlotInfo(fig, plotInfo)
        fig.canvas.draw()

        return fig


def addPlotInfo(fig: Figure, plotInfo: Mapping[str, str]) -> Figure:
    """Add useful information to the plot.

    Parameters
    ----------
    fig : `matplotlib.figure.Figure`
        The figure to add the information to.
    plotInfo : `dict`
        A dictionary of the plot information.

    Returns
    -------
    fig : `matplotlib.figure.Figure`
        The figure with the information added.
    """
    fig.text(0.01, 0.99, plotInfo["plotName"], fontsize=7, transform=fig.transFigure, ha="left", va="top")
    infoText = parsePlotInfo(plotInfo)
    fig.text(0.01, 0.984, infoText, fontsize=6, transform=fig.transFigure, alpha=0.6, ha="left", va="top")

    return fig


def parsePlotInfo(plotInfo: Mapping[str, str]) -> str:
    """Extract information from the plotInfo dictionary and parses it into
    a meaningful string that can be added to a figure. The default function
    in .plotUtils is not suitable for image plotting.

    Parameters
    ----------
    plotInfo : `dict`[`str`, `str`]
        A plotInfo dictionary containing useful information to
        be included on a figure.

    Returns
    -------
    infoText : `str`
        A string containing the plotInfo information, parsed in such a
        way that it can be included on a figure.
    """
    run = plotInfo["run"]
    componentType = f"\nComponent: {plotInfo['component']}"

    maskPlaneText = ""
    if "maskPlanes" in plotInfo:
        for maskPlane in plotInfo["maskPlanes"]:
            maskPlaneText += maskPlane + ", "
            maskPlaneText = f", Mask Plane(s): {maskPlaneText[:-2]}"

    dataIdText = f"\nSkyMap:{plotInfo['skymap']}, Tract: {plotInfo['tract']}"

    bandText = ""
    for band in plotInfo["bands"]:
        bandText += band + ", "
    bandsText = f", Bands: {bandText[:-2]}"
    infoText = f"\n{run}{componentType}{maskPlaneText}{dataIdText}{bandsText}"

    return infoText
