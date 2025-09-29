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
    "PerTractPropertyMapPlot",
    "SurveyWidePropertyMapPlot",
)

import logging
from typing import Mapping, Union

import lsst.pex.config as pexConfig
import lsst.sphgeom as sphgeom
import matplotlib.patheffects as mpl_path_effects
import numpy as np
import skyproj
from healsparse.healSparseMap import HealSparseMap
from lsst.analysis.tools.tasks.propertyMapAnalysis import (
    PerTractPropertyMapAnalysisConfig,
    SurveyWidePropertyMapAnalysisConfig,
)
from lsst.skymap.tractInfo import ExplicitTractInfo
from lsst.utils.plotting import make_figure, set_rubin_plotstyle
from matplotlib import cm, rc_context
from matplotlib.figure import Figure
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_size import AxesX, AxesY, Fraction

from ...interfaces import KeyedData, PlotAction

_LOG = logging.getLogger(__name__)

# Holds unit renames to match style guidelines: {"old_unit": "new_unit"}.
unitRenameDict = {"mag(AB)": r"$\rm mag_{AB}$"}


def getZoomedExtent(box, n):
    """Get a zoomed plot extent in degrees for a bounding box.

    Parameters
    ----------
    box : `lsst.sphgeom.NormalizedAngleInterval`
        The bounding box to get an extent for.
    n : `float`
        Zoom factor; for instance, n=2 means zooming in 2 times at the
        center. Must be positive.

    Returns
    -------
    `tuple` [`float`]
        New extent as (new_lon_min, new_lon_max, new_lat_min, new_lat_max).
    """
    lon = box.getLon()
    lat = box.getLat()
    if n > 1:
        factor = (1.0 - 1.0 / n) / 2.0
        lon.erodeBy(sphgeom.NormalizedAngle(factor * box.getWidth()))
        lat.erodeBy(sphgeom.NormalizedAngle(factor * box.getHeight()))
    elif n < 1:
        factor = (1.0 / n - 1) / 2.0
        lon.dilateBy(sphgeom.NormalizedAngle(factor * box.getWidth()))
        lat.dilateBy(sphgeom.NormalizedAngle(factor * box.getHeight()))

    extent = (
        lon.getA().asDegrees(),
        lon.getB().asDegrees(),
        lat.getA().asDegrees(),
        lat.getB().asDegrees(),
    )
    return extent


def getLongestSuffixMatch(s, options):
    """Find the longest suffix in the provided list that matches the end of
    the given string.

    Parameters
    ----------
    s : `str`
        The target string for which we want to find a matching suffix.
    options : `list` [`str`]
        A list of potential suffix strings to match against the target
        string `s`.

    Returns
    -------
    `str`
        The longest matching suffix from the `options` list. If no match is
        found, returns `None`.
    """
    return next((opt for opt in sorted(options, key=len, reverse=True) if s.endswith(opt)), None)


def addTextToColorbar(
    cb, text, orientation="vertical", color="black", fontsize=14, fontweight="bold", alpha=0.8
):
    """Helper method to add text inside the horizontal colorbar.

    Parameters
    ----------
    cb : `~matplotlib.colorbar.Colorbar`
        The colorbar object.
    text : `str`
        The text to add.
    orientation : `str`, optional
        The orientation of the colorbar. Can be either "vertical" or
        "horizontal".
    fontsize : `int`, optional
        The fontsize of the text.
    fontweight : `str`, optional
        The fontweight of the text.
    alpha : `float`, optional
        The alpha value of the text.

    Returns
    -------
    `None`
        The text is added to the colorbar in place.
    """
    if color is None:
        color = "black"
    vmid = (cb.vmin + cb.vmax) / 2
    positions = {"vertical": (0.5, vmid), "horizontal": (vmid, 0.5)}
    cbtext = cb.ax.text(
        *positions[orientation],
        text,
        color=color,
        va="center",
        ha="center",
        fontsize=fontsize,
        fontweight=fontweight,
        rotation=orientation,
        alpha=alpha,
    )
    # Add a distinct outline around the text for better visibility in
    # various backgrounds.
    cbtext.set_path_effects(
        [mpl_path_effects.Stroke(linewidth=4, foreground="white", alpha=0.8), mpl_path_effects.Normal()]
    )


class CustomHandler(HandlerTuple):
    """Custom legend handler to overlay multiple patches for a single
    legend entry.

    This handler class inherits from `HandlerTuple` and is designed to
    handle cases where multiple artists (e.g., patches) need to be overlaid
    on top of each other in a single legend entry, as opposed to
    side-by-side which is the default behavior of `HandlerTuple`.

    Methods
    -------
    create_artists:
        Override the `create_artists` method of `HandlerTuple` to modify
        the positioning of the artists so that they overlay directly on top
        of one another in the legend.

    Example
    -------
    # Plot some data.
    line, = ax.plot(x, y, label="Sample Line")

    # Use CustomHandler for overlaid patches and also include the regular
    # line legend if desired.
    handles = [(patch1, patch2), line]
    labels = ['Overlaid Patches', line.get_label()]
    leg = ax.legend(
        handles, labels, handler_map={tuple: CustomHandler()}, loc="best"
    )
    """

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        artists = HandlerTuple.create_artists(
            self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
        )
        # Overlay the two patches.
        for a in artists:
            a.set_transform(trans)
        return artists


class PerTractPropertyMapPlot(PlotAction):
    draw_patch_bounds = pexConfig.Field[bool](
        doc="Whether to draw patch inner boundaries or not",
        default=False,
    )
    label_patches = pexConfig.ChoiceField[str](
        doc="Which patches to label by ID",
        default="none",
        allowed={
            "all": "All patches",
            "edge": "Edge patches only",
            "none": "No labels",
        },
    )
    plotName = pexConfig.Field[str](doc="The name for the plotting task.", optional=True)

    def __call__(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        plotConfig: PerTractPropertyMapAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
        **kwargs,
    ) -> Mapping[str, Figure]:
        self._validateInput(data, tractInfo, plotConfig, plotInfo)
        return self.makePlot(data, tractInfo, plotConfig, plotInfo)

    def _validateInput(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        plotConfig: PerTractPropertyMapAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
    ) -> None:
        """Validate the input data."""

        if not isinstance(tractInfo, ExplicitTractInfo):
            raise TypeError(f"Input `tractInfo` type must be {ExplicitTractInfo} not {type(tractInfo)}.")

        if not isinstance(plotConfig, PerTractPropertyMapAnalysisConfig):
            raise TypeError(
                "`plotConfig` must be a "
                "`lsst.analysis.tools.tasks.propertyMapAnalysis.PerTractPropertyMapAnalysisConfig`. "
                f"Got {type(plotConfig)}."
            )

        if not isinstance(plotInfo, dict):
            raise TypeError("`plotConfig` must be a dictionary.")

        if not plotConfig.publicationStyle:
            zoomFactors = plotConfig.zoomFactors
            isListOfFloats = isinstance(zoomFactors, pexConfig.listField.List) and all(
                isinstance(zf, float) for zf in zoomFactors
            )
            if not (isListOfFloats and len(zoomFactors) == 2) or any(zf <= 1 for zf in zoomFactors):
                raise TypeError(
                    "`zoomFactors` must be a two-element `lsst.pex.config.listField.List` of floats > 1."
                )

        for atool in plotConfig.atools:
            if not isinstance(atool.nBinsHist, int) or atool.nBinsHist <= 0:
                raise ValueError(
                    f"`nBinsHist` for property `{atool.process.buildActions.data.mapKey}` must be a positive "
                    f"integer. Got {atool.nBinsHist}."
                )

        # Identify any invalid entries in `data`.
        invalidEntries = {
            key: pytype
            for key, pytype in {k: v.ref.datasetType.storageClass.pytype for k, v in data.items()}.items()
            if pytype != HealSparseMap
        }

        # If any invalid entries are found, raise a TypeError with details.
        if invalidEntries:
            errorMessage = "; ".join(
                f"`{key}` should be {HealSparseMap}, got {type_}" for key, type_ in invalidEntries.items()
            )
            raise TypeError(f"Invalid input types found in `data`: {errorMessage}")

    def addPlotInfo(
        self,
        fig: Figure,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
        toolName: str,
    ) -> Figure:
        """Add useful information to the plot.

        Parameters
        ----------
        fig : `matplotlib.figure.Figure`
            The figure to add the information to.
        plotInfo : `dict`
            A dictionary of the plot information.
        toolName : `str`
            The name of the tool used to generate the plot.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The figure with the information added.
        """

        run = plotInfo["run"]
        tableType = f"\nTable: {plotInfo['tableNames'][toolName]}"

        dataIdText = f"Tract: {plotInfo['tract']}, Band: {plotInfo['band']}"
        propertyDescription = plotInfo["description"]
        # Lowercase the first letter unless the string starts with more than
        # one uppercase letter in which case we keep it as is, e.g. PSF, DCR.
        if propertyDescription[0].isupper() and propertyDescription[1].islower():
            propertyDescription = propertyDescription[0].lower() + propertyDescription[1:]
        mapText = (
            f", Property: {propertyDescription}, "
            f"Unit: {plotInfo['unit']}, "
            f"Operation: {plotInfo['operation']}, "
            f"Coadd: {plotInfo['coaddName']}"
        )
        geomText = (
            f", Valid area: {plotInfo['valid_area']:.2f} sq. deg., "
            + f"NSIDE: {plotInfo['nside']}, projection: {plotInfo['projection']}"
        )
        infoText = f"\n{dataIdText}{mapText}"

        fig.text(
            0.04,
            0.965,
            f'{plotInfo["plotName"]} of {plotInfo["property"]}',
            fontsize=19,
            transform=fig.transFigure,
            ha="left",
            va="top",
        )
        t = fig.text(
            0.04,
            0.942,
            f"{run}{tableType}{geomText}{infoText}",
            fontsize=15,
            transform=fig.transFigure,
            alpha=0.6,
            ha="left",
            va="top",
        )
        t.set_linespacing(1.4)

        return fig

    def makePlot(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        plotConfig: PerTractPropertyMapAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
    ) -> Mapping[str, Figure]:
        """Make the survey property map plot.

        Parameters
        ----------
        data : `KeyedData`
            The HealSparseMap to plot the points from.
        tractInfo: `~lsst.skymap.tractInfo.ExplicitTractInfo`
            The tract info object.
        plotConfig :
            `~lsst.analysis.tools.tasks.perTractPropertyMapAnalysis.
            PerTractPropertyMapAnalysisConfig`
            The configuration for the plot.
        plotInfo : `dict`
            A dictionary of information about the data being plotted.

        Returns
        -------
        fig : `~matplotlib.figure.Figure`
            The resulting figure.
        """

        set_rubin_plotstyle()

        # 'plotName' by default is constructed from the attribute specified in
        # 'atools.<attribute>' in the pipeline YAML. If it is explicitly
        # set in `~lsst.analysis.tools.atools.propertyMap.PropertyMapTool`,
        # it will override this default.
        if self.plotName:
            # Set the plot name using 'produce.plot.plotName' from
            # PropertyMapTool's instance.
            plotInfo["plotName"] = self.plotName

        # Plotting customization.
        colorbarLocation = "right"
        colorbarTickLabelSize = 14
        colorBarAspect = 16
        colorBarFraction = 0.07  # w.r.t. the axes size.
        colorBarPad = 0  # w.r.t. the axes size.
        rcparams = {
            "axes.labelsize": 22 if plotConfig.publicationStyle else 19,
            "axes.linewidth": 1.8,
            "xtick.labelsize": 15 if plotConfig.publicationStyle else 14,
            "ytick.labelsize": 15 if plotConfig.publicationStyle else 14,
        }

        colorbarKwargs = dict(plotConfig.colorbarKwargs)
        projectionKwargs = dict(plotConfig.projectionKwargs)

        if plotConfig.publicationStyle:
            zoomFactors = []
            _LOG.info(
                "Zoom factors are not used in publication-style plots. "
                "Only the full-tract map will be plotted."
            )
            if "cmap" in colorbarKwargs and colorbarKwargs["cmap"] != "viridis":
                _LOG.warning(
                    "Color map set to 'viridis' for publication style plots. "
                    f"The color map '{colorbarKwargs['cmap']}' set in the config will be ignored."
                )
            colorbarKwargs["cmap"] = "viridis"
            cmap = cm.get_cmap(colorbarKwargs["cmap"])
            # Colorbar text color only, not used for any histogram.
            histColors = [cmap(0.1)]
            labelpad = 14
        else:
            zoomFactors = plotConfig.zoomFactors
            colorbarKwargs["cmap"] = colorbarKwargs.get("cmap", "viridis")
            # Muted green for the full map, and muted red and blue for the two
            # zoomed-in maps. Used for boxes, colorbar texts and histograms.
            histColors = ["#265D40", "#8B0000", "#00008B"]
            labelpad = 11

        toolName = data["data"].ref.datasetType.name
        mapName = toolName.replace("_map_", "_")
        mapData = data["data"].get()

        with rc_context(rcparams):
            if plotConfig.publicationStyle:
                # Make a single plot for the full tract.
                fig = make_figure(figsize=(8, 8))
                ax1 = fig.add_subplot(111)
                fig.subplots_adjust(left=0.133, right=0.888, bottom=0.1, top=0.855)
            else:
                # Make a 2x2 grid of subplots for the full tract, two zoomed-in
                # views, and a histogram of values.
                fig = make_figure(figsize=(16, 16))
                ax1 = fig.add_subplot(221)
                ax2 = fig.add_subplot(222)
                ax3 = fig.add_subplot(223)
                ax4 = fig.add_subplot(224)

                # Reduce whitespace but leave some room at the top for info.
                fig.subplots_adjust(left=0.064, right=0.96, top=0.855, bottom=0.07, wspace=0.275, hspace=0.24)

            # Get the values for the valid pixels of the full tract.
            values = mapData[mapData.valid_pixels]
            goodValues = np.isfinite(values)
            values = values[goodValues]  # As a precaution.

            # Make a concise human-readable label for the plot.
            plotInfo["unit"] = "N/A"  # Unless overridden.
            if hasattr(mapData, "metadata") and all(
                key in mapData.metadata for key in ["DESCRIPTION", "OPERATION", "UNIT"]
            ):
                hasMetadata = True
                metadata = mapData.metadata
                plotInfo["description"] = metadata["DESCRIPTION"]
                plotInfo["operation"] = metadata["OPERATION"]
                if metadata["UNIT"]:
                    plotInfo["unit"] = metadata["UNIT"]
                elif metadata["UNIT"] == "":
                    plotInfo["unit"] = "dimensionless"
            else:
                hasMetadata = False
                plotInfo["operation"] = getLongestSuffixMatch(
                    mapName, ["min", "max", "mean", "weighted_mean", "sum"]
                ).replace("_", " ")
            plotInfo["coaddName"] = mapName.split("Coadd_")[0]
            plotInfo["operation"] = plotInfo["operation"].replace("minimum", "min").replace("maximum", "max")
            propertyName = mapName[len(f"{plotInfo['coaddName']}Coadd_") : -len(plotInfo["operation"])].strip(
                "_"
            )
            if not hasMetadata:
                # Infer the property description from the map name (all
                # lower case), and properly handle formatting.
                plotInfo["description"] = (
                    propertyName.replace("_", " ")
                    .replace("psf", "PSF")
                    .replace("dcr", "DCR")
                    .replace("dra", r"$\Delta$RA")
                    .replace("ddec", r"$\Delta$Dec")
                )
            plotInfo["property"] = (
                propertyName.replace("_", " ")
                .title()  # Capitalize and handle edge cases below.
                .replace("Psf", "PSF")
                .replace("Dcr", "DCR")
                .replace("Dra", r"$\Delta$RA")
                .replace("Ddec", r"$\Delta$Dec")
                .replace("E1", "e1")
                .replace("E2", "e2")
            )

            # Handle unit renaming.
            if plotInfo["unit"] in unitRenameDict:
                plotInfo["unit"] = unitRenameDict[plotInfo["unit"]]

            atool = getattr(plotConfig.atools, toolName)
            nBinsHist = atool.nBinsHist

            # Use the full tract bounding box to set the default extent.
            box = tractInfo.getOuterSkyPolygon().getBoundingBox()
            ctr_lonlat = box.getCenter()
            ctr_lon = ctr_lonlat.getLon().asDegrees()
            ctr_lat = ctr_lonlat.getLat().asDegrees()

            # Prepare plot elements for the full tract and optional zoomed-in
            # views.
            if plotConfig.publicationStyle:
                skyprojAxes = [ax1]
            else:
                skyprojAxes = [ax1, ax3, ax4]

            zoomIdx = []
            for ax, zoomFactor, histColor in zip(skyprojAxes, [1.0, *zoomFactors], histColors):
                extent = getZoomedExtent(box, zoomFactor)

                sp = skyproj.GnomonicSkyproj(
                    ax=ax,
                    lon_0=ctr_lon,
                    lat_0=ctr_lat,
                    extent=extent,
                    **projectionKwargs,
                )

                sp.draw_hspmap(mapData, zoom=False, cmap=colorbarKwargs["cmap"])

                sp.ax.set_xlabel("R.A.", labelpad=labelpad, fontsize=rcparams["axes.labelsize"])
                sp.ax.set_ylabel("Dec.", labelpad=labelpad, fontsize=rcparams["axes.labelsize"])

                if self.draw_patch_bounds:
                    label_all = self.label_patches == "all"
                    label_edge = self.label_patches == "edge"
                    if label_all or label_edge:
                        patchids_max = tuple(x - 1 for x in tractInfo.getNumPatches())
                        lon_min, lon_max = min(extent[:2]), max(extent[:2])
                        lat_min, lat_max = min(extent[2:]), max(extent[2:])

                    for patchInfo in tractInfo:
                        vertices = patchInfo.getInnerSkyPolygon().getVertices()
                        clipped = tractInfo.inner_sky_region.clipTo(
                            sphgeom.Box(sphgeom.LonLat(vertices[0]), sphgeom.LonLat(vertices[2]))
                        )
                        lonlats_patch = np.array(
                            [
                                [x.asDegrees() for x in (lonlat.getA(), lonlat.getB())]
                                for lonlat in (clipped.getLon(), clipped.getLat())
                            ]
                        )
                        lons_patch = np.concat((lonlats_patch[0, :], lonlats_patch[0, ::-1]))
                        lats_patch = np.repeat(lonlats_patch[1, :], 2)

                        sp.draw_polygon(lons_patch, lats_patch, edgecolor="gray")
                        label_id = label_all or (
                            label_edge
                            and (
                                (patchInfo.index[0] == 0)
                                or (patchInfo.index[1] == 0)
                                or (patchInfo.index[0] == patchids_max[0])
                                or (patchInfo.index[1] == patchids_max[1])
                            )
                        )
                        if label_id:
                            lon_label = clipped.getCenter().getLon().asDegrees()
                            lat_label = clipped.getCenter().getLat().asDegrees()
                            if (lon_min < lon_label < lon_max) and (lat_min < lat_label < lat_max):
                                sp.ax.text(
                                    lon_label,
                                    lat_label,
                                    patchInfo.getSequentialIndex(),
                                    ha="center",
                                    va="center",
                                    size=7,
                                    c=[0, 0, 0, 0.5],
                                )

                # Specify the size and padding of the colorbar axes with
                # respect to the main axes.
                refAx = AxesX(ax) if colorbarLocation in ("left", "right") else AxesY(ax)
                cbsize = Fraction(colorBarFraction, refAx)
                cbpad = Fraction(colorBarPad, refAx)

                # Make divider and a colorbar axes to be attached to the main
                # axes.
                divider = make_axes_locatable(ax)
                cax = divider.append_axes(colorbarLocation, size=cbsize, pad=cbpad)

                cbar = sp.draw_colorbar(
                    **{
                        "location": colorbarLocation,
                        "aspect": colorBarAspect,
                        "fraction": colorBarFraction,
                        "pad": colorBarPad,
                        "cax": cax,
                        **colorbarKwargs,
                    }
                )
                cbar.ax.tick_params(labelsize=colorbarTickLabelSize)
                if plotConfig.publicationStyle:
                    if plotInfo["unit"] not in ["dimensionless", "N/A"]:
                        cbarText = f"{plotInfo['property']} ({plotInfo['unit']})"
                    else:
                        cbarText = plotInfo["property"]
                else:
                    cbarText = (
                        "Full Tract" if zoomFactor == 1.0 else f"{self.prettyPrintFloat(zoomFactor)}x Zoom"
                    )
                addTextToColorbar(cbar, cbarText, color=histColor)
                if zoomFactor == 1.0:
                    # Store the "full tract" map so that we can overplot
                    # the zoom rectangles.
                    spf = sp
                else:
                    # Create a rectangle for the zoomed-in region.
                    x0, x1, y0, y1 = extent
                    for c, ls, lw in zip(["white", histColor], ["solid", "dashed"], [3.5, 1.5]):
                        spf.draw_polygon(
                            [x0, x0, x1, x1],
                            [y0, y1, y1, y0],
                            facecolor="none",
                            edgecolor=c,
                            linestyle=ls,
                            linewidth=lw,
                            alpha=0.8,
                        )
                    zoomText = spf.ax.text(
                        (x0 + x1) / 2,
                        y0,
                        f"{self.prettyPrintFloat(zoomFactor)}x",
                        color=histColor,
                        fontsize=14,
                        fontweight="bold",
                        alpha=0.8,
                        ha="center",
                        va="bottom",
                    )
                    # Add a distinct outline around the text for better
                    # visibility in various backgrounds.
                    zoomText.set_path_effects(
                        [
                            mpl_path_effects.Stroke(linewidth=4, foreground="white", alpha=0.8),
                            mpl_path_effects.Normal(),
                        ]
                    )
                    # Get the indices of pixels in the zoomed-in region.
                    pos = mapData.valid_pixels_pos()
                    # Reversed axes consideration.
                    xmin, xmax = sorted([x0, x1])
                    idx = (pos[0] > xmin) & (pos[0] < xmax) & (pos[1] > y0) & (pos[1] < y1)
                    zoomIdx.append(idx[goodValues])

            # Calculate weights for each bin to ensure that the peak of the
            # histogram reaches 1.
            if len(values) == 0 or np.ptp(values) < 1e-12:
                # np.histogram cannot use 100 bins if no variance in values
                weights = np.ones_like(values) / values.size
                nBinsHist = 1
                _LOG.info(f"No variance in {toolName}; set histogram to a single bin.")
            else:
                weights = np.ones_like(values) / np.histogram(values, bins=nBinsHist)[0].max()

            if not plotConfig.publicationStyle:
                # Compute full-tract histogram and get its bins.
                # NOTE: `exposure_time` histograms are quantized and look more
                # bar-like, so they look better with fewer bins.
                bins = ax2.hist(
                    values,
                    bins=nBinsHist,
                    label="Full Tract",
                    color=histColors[0],
                    weights=weights,
                    alpha=0.7,
                )[1]

                # Align the histogram (top right panel) with the skyproj plots.
                pos1 = spf.ax.get_position()  # Top left.
                pos4 = sp.ax.get_position()  # Bottom right.
                cbarWidth = cbar.ax.get_position().height / colorBarAspect
                # NOTE: cbarWidth != cbar.ax.get_position().width
                ax2.set_position([pos4.x0, pos1.y0, pos4.width + cbarWidth, pos1.height])

                # Overplot the histograms for the zoomed-in plots.
                for zoomFactor, zidx, color, linestyle, hatch in zip(
                    zoomFactors, zoomIdx, histColors[1:], ["solid", "dotted"], ["//", "xxxx"]
                ):
                    weights = np.ones_like(values[zidx]) / np.histogram(values[zidx], bins=bins)[0].max()
                    histLabel = f"{self.prettyPrintFloat(zoomFactor)}x Zoom"
                    histValues = ax2.hist(
                        values[zidx],
                        bins=bins,
                        label=histLabel,
                        color=color,
                        weights=weights,
                        histtype="step",
                        linewidth=2,
                        linestyle=linestyle,
                        alpha=0.6,
                    )[0]
                    # Fill the area under the step.
                    ax2.fill_between(
                        bins[:-1],
                        histValues,
                        step="post",
                        color=color,
                        alpha=0.2,
                        hatch=hatch,
                        label="hidden",
                    )

                # Set labels and legend.
                xlabel = plotInfo["property"]
                if plotInfo["unit"] not in ["dimensionless", "N/A"]:
                    xlabel += f" ({plotInfo['unit']})"
                xtext = ax2.set_xlabel(xlabel, labelpad=labelpad)
                ytext = ax2.set_ylabel("Normalized Count", labelpad=labelpad)
                xtext.set_fontsize(rcparams["axes.labelsize"])
                ytext.set_fontsize(rcparams["axes.labelsize"])

                # Get handles and labels from the axis.
                handles, labels = ax2.get_legend_handles_labels()

                # Add a legend with custom handler that combines the handle
                # pairs for the zoomed-in cases.
                handles = [handles[0], (handles[1], handles[2]), (handles[3], handles[4])]
                while "hidden" in labels:
                    labels.remove("hidden")
                legend = ax2.legend(
                    handles,
                    labels,
                    handler_map={tuple: CustomHandler()},
                    loc="best",
                    frameon=False,
                    fontsize=15,
                )

                for line, text in zip(handles, legend.get_texts()):
                    if isinstance(line, tuple):
                        # Use the first handle to get the color.
                        line = line[0]
                    color = line.get_edgecolor() if line.get_facecolor()[-1] == 0 else line.get_facecolor()
                    text.set_color(color)

                # Add extra info to plotInfo.
                plotInfo["projection"] = "Gnomonic"
                plotInfo["nside"] = mapData.nside_sparse
                plotInfo["valid_area"] = mapData.get_valid_area()

                # Add useful information to the plot.
                self.addPlotInfo(fig, plotInfo, toolName)
                style = ""
            else:
                style = "publication-style "
                fig.suptitle(
                    f"{plotInfo['description']} {plotInfo['operation']} map",
                    fontsize=18.5,
                    ha="center",
                    va="top",
                    y=0.985,
                )
                fig.text(
                    0.5,
                    0.925,
                    f"Tract: {plotInfo['tract']}, Band: {plotInfo['band']}, Coadd: {plotInfo['coaddName']}",
                    ha="center",
                    fontsize=15.5,
                )

            _LOG.info(
                f"Made {style}per-tract property map plot for dataset type '{toolName}', "
                f"tract: {plotInfo['tract']}, band: '{plotInfo['band']}'."
            )

        return fig

    @staticmethod
    def prettyPrintFloat(n):
        if n.is_integer():
            return str(int(n))
        return str(n)


class SurveyWidePropertyMapPlot(PlotAction):
    plotName = pexConfig.Field[str](doc="The name for the plotting task.", optional=True)

    def __call__(
        self,
        data: KeyedData,
        plotConfig: SurveyWidePropertyMapAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
        **kwargs,
    ) -> Mapping[str, Figure]:
        return self.makePlot(data, plotConfig, plotInfo)

    def addPlotInfo(
        self,
        fig: Figure,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
        toolName: str,
    ) -> Figure:
        """Add useful information to the plot.

        Parameters
        ----------
        fig : `matplotlib.figure.Figure`
            The figure to add the information to.
        plotInfo : `dict`
            A dictionary of the plot information.
        toolName : `str`
            The name of the tool used to generate the plot.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The figure with the information added.
        """

        run = plotInfo["run"]
        tableType = f"\nTable: {plotInfo['tableNames'][toolName]}"

        dataIdText = f"Band: {plotInfo['band']}"
        propertyDescription = plotInfo["description"]
        # Lowercase the first letter unless the string starts with more than
        # one uppercase letter in which case we keep it as is, e.g. PSF, DCR.
        if propertyDescription[0].isupper() and propertyDescription[1].islower():
            propertyDescription = propertyDescription[0].lower() + propertyDescription[1:]
        mapText = (
            f", Property: {propertyDescription}, "
            f"Unit: {plotInfo['unit']}, "
            f"Operation: {plotInfo['operation']}, "
            f"Coadd: {plotInfo['coaddName']}"
        )
        geomText = (
            f", Valid area: {plotInfo['valid_area']:.2f} sq. deg., "
            + f"NSIDE: {plotInfo['nside']}, projection: {plotInfo['projection']}"
        )
        infoText = f"\n{dataIdText}{mapText}"

        titleBoxTopLeftCorner = (0.045, 0.89)
        title = fig.text(
            *titleBoxTopLeftCorner,
            f'{plotInfo["plotName"]} of {plotInfo["property"]}',
            fontsize=19,
            transform=fig.transFigure,
            ha="left",
            va="top",
        )
        lineHeightFraction = title.get_fontsize() / (fig.get_size_inches()[1] * fig.dpi)
        infoBoxTopLeftCorner = (titleBoxTopLeftCorner[0], titleBoxTopLeftCorner[1] - 1.8 * lineHeightFraction)
        info = fig.text(
            *infoBoxTopLeftCorner,
            f"{run}{tableType}{geomText}{infoText}",
            fontsize=15,
            transform=fig.transFigure,
            alpha=0.6,
            ha="left",
            va="top",
        )
        info.set_linespacing(1.4)

        return fig

    def makePlot(
        self,
        data: KeyedData,
        plotConfig: SurveyWidePropertyMapAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
    ) -> Figure:
        """Make the survey property map plot.

        Parameters
        ----------
        data : `KeyedData`
            The HealSparseMap to plot the points from.
        plotConfig :
            `~lsst.analysis.tools.tasks.propertyMapSurveyAnalysis.
            SurveyWidePropertyMapAnalysisConfig`
            The configuration for the plot.
        plotInfo : `dict`
            A dictionary of information about the data being plotted.

        Returns
        -------
        fig : `~matplotlib.figure.Figure`
            The resulting figure.
        """

        set_rubin_plotstyle()

        # 'plotName' by default is constructed from the attribute specified in
        # 'atools.<attribute>' in the pipeline YAML. If it is explicitly
        # set in `~lsst.analysis.tools.atools.healSparsePropertyMap.
        # SurveyWidePropertyMapTool`, it will override this default.
        if self.plotName:
            # Set the plot name using 'produce.plot.plotName' from
            # SurveyWidePropertyMapTool's instance.
            plotInfo["plotName"] = self.plotName

        # Plotting customization.
        colorbarTickLabelSize = 16
        rcparams = {
            "axes.labelsize": 18,
            "axes.linewidth": 1.8,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
        }

        toolName = data["data"].ref.datasetType.name
        mapName = toolName.replace("_consolidated_map_", "_")
        mapData = data["data"].get()

        with rc_context(rcparams):
            # The figsize should be decided based on the survey footprint.
            fig = make_figure(figsize=(19, 7))
            ax = fig.add_subplot(111)

            if not plotConfig.publicationStyle:
                # Leave some room at the top for plotInfo.
                fig.subplots_adjust(left=0.072, right=0.945, top=0.55)

            # Get the values for the valid pixels of the full tract.
            values = mapData[mapData.valid_pixels]
            goodValues = np.isfinite(values)
            values = values[goodValues]  # As a precaution.

            # Make a concise human-readable label for the plot.
            plotInfo["unit"] = "N/A"  # Unless overridden.
            if hasattr(mapData, "metadata") and all(
                key in mapData.metadata for key in ["DESCRIPTION", "OPERATION", "UNIT"]
            ):
                hasMetadata = True
                metadata = mapData.metadata
                plotInfo["description"] = metadata["DESCRIPTION"]
                plotInfo["operation"] = metadata["OPERATION"]
                if metadata["UNIT"]:
                    plotInfo["unit"] = metadata["UNIT"]
                elif metadata["UNIT"] == "":
                    plotInfo["unit"] = "dimensionless"
            else:
                hasMetadata = False
                plotInfo["operation"] = getLongestSuffixMatch(
                    mapName, ["min", "max", "mean", "weighted_mean", "sum"]
                ).replace("_", " ")
            plotInfo["coaddName"] = mapName.split("Coadd_")[0]
            plotInfo["operation"] = plotInfo["operation"].replace("minimum", "min").replace("maximum", "max")
            propertyName = mapName[
                len(f"{plotInfo['coaddName']}Coadd_") : -len(f"{plotInfo['operation']}")
            ].strip("_")
            if not hasMetadata:
                # Infer the property description from the map name (all
                # lower case), and properly handle formatting.
                plotInfo["description"] = (
                    propertyName.replace("_", " ")
                    .replace("psf", "PSF")
                    .replace("dcr", "DCR")
                    .replace("dra", r"$\Delta$RA")
                    .replace("ddec", r"$\Delta$Dec")
                )
            plotInfo["property"] = (
                propertyName.replace("_", " ")
                .title()  # Capitalize and handle edge cases below.
                .replace("Psf", "PSF")
                .replace("Dcr", "DCR")
                .replace("Dra", r"$\Delta$RA")
                .replace("Ddec", r"$\Delta$Dec")
                .replace("E1", "e1")
                .replace("E2", "e2")
            )

            # Handle unit renaming.
            if plotInfo["unit"] in unitRenameDict:
                plotInfo["unit"] = unitRenameDict[plotInfo["unit"]]

            sp = getattr(skyproj, f"{plotConfig.projection}Skyproj")(ax=ax, **plotConfig.projectionKwargs)

            colorbarKwargs = dict(plotConfig.colorbarKwargs)
            if plotConfig.publicationStyle:
                colorbarKwargs["cmap"] = "viridis"
            else:
                colorbarKwargs["cmap"] = colorbarKwargs.get("cmap", "viridis")

            # Work around skyproj bug that will fail to zoom on empty map.
            if mapData.n_valid == 0:
                if plotConfig.autozoom:
                    _LOG.warning("No valid pixels found in the map. Auto zooming is disabled.")
                sp.draw_hspmap(mapData, zoom=False, cmap=colorbarKwargs["cmap"])
            else:
                sp.draw_hspmap(mapData, zoom=plotConfig.autozoom, cmap=colorbarKwargs["cmap"])
            sp.ax.set_xlabel("R.A.")
            sp.ax.set_ylabel("Dec.")

            # In the below, colorbarKwargs takes precedence over hardcoded
            # arguments in case of conflict.
            cbar = sp.draw_colorbar(**{"location": "top", "pad": 0.2, **colorbarKwargs})
            cbar.ax.tick_params(labelsize=colorbarTickLabelSize)
            unit = f" ({plotInfo['unit']})" if plotInfo["unit"] not in ["dimensionless", "N/A"] else ""
            cbarText = f"{plotInfo['property']}{unit}"
            cbarLoc = colorbarKwargs["location"]
            cbarOrientation = colorbarKwargs.get("orientation", None)
            if cbarOrientation is None:
                cbarOrientation = "vertical" if cbarLoc in ["right", "left"] else "horizontal"
            addTextToColorbar(cbar, cbarText, color="#265D40", fontsize=16, orientation=cbarOrientation)

            # Add extra info to plotInfo.
            plotInfo["projection"] = plotConfig.projection
            plotInfo["nside"] = mapData.nside_sparse
            plotInfo["valid_area"] = mapData.get_valid_area()

            # Add useful information to the plot.
            self.addPlotInfo(fig, plotInfo, toolName)

        _LOG.info(
            f"Made survey-wide property map plot for dataset type '{toolName}', "
            f"band: '{plotInfo['band']}'."
        )

        return fig
