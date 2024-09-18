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

__all__ = ("PanelConfig",)

from typing import TYPE_CHECKING, Iterable, List, Mapping, Tuple

import matplotlib
import numpy as np
from lsst.geom import Box2D, SpherePoint, degrees
from lsst.pex.config import Config, Field
from matplotlib import cm, colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from scipy.stats import binned_statistic_2d

from ...math import nanMedian, nanSigmaMad

if TYPE_CHECKING:
    from matplotlib.figure import Figure

null_formatter = matplotlib.ticker.NullFormatter()


def generateSummaryStats(data, skymap, plotInfo):
    """Generate a summary statistic in each patch or detector.

    Parameters
    ----------
    data : `dict`
        A dictionary of the data to be plotted.
    skymap : `lsst.skymap.BaseSkyMap`
        The skymap associated with the data.
    plotInfo : `dict`
        A dictionary of the plot information.

    Returns
    -------
    patchInfoDict : `dict`
        A dictionary of the patch information.
    """
    tractInfo = skymap.generateTract(plotInfo["tract"])
    tractWcs = tractInfo.getWcs()

    # For now also convert the gen 2 patchIds to gen 3
    if "y" in data.keys():
        yCol = "y"
    elif "yStars" in data.keys():
        yCol = "yStars"
    elif "yGalaxies" in data.keys():
        yCol = "yGalaxies"
    elif "yUnknowns" in data.keys():
        yCol = "yUnknowns"

    patchInfoDict = {}
    maxPatchNum = tractInfo.num_patches.x * tractInfo.num_patches.y
    patches = np.arange(0, maxPatchNum, 1)
    for patch in patches:
        if patch is None:
            continue
        # Once the objectTable_tract catalogues are using gen 3 patches
        # this will go away
        onPatch = data["patch"] == patch
        if sum(onPatch) == 0:
            stat = np.nan
        else:
            stat = nanMedian(data[yCol][onPatch])
        try:
            patchTuple = (int(patch.split(",")[0]), int(patch.split(",")[-1]))
            patchInfo = tractInfo.getPatchInfo(patchTuple)
            gen3PatchId = tractInfo.getSequentialPatchIndex(patchInfo)
        except AttributeError:
            # For native gen 3 tables the patches don't need converting
            # When we are no longer looking at the gen 2 -> gen 3
            # converted repos we can tidy this up
            gen3PatchId = patch
            patchInfo = tractInfo.getPatchInfo(patch)

        corners = Box2D(patchInfo.getInnerBBox()).getCorners()
        skyCoords = tractWcs.pixelToSky(corners)

        patchInfoDict[gen3PatchId] = (skyCoords, stat)

    tractCorners = Box2D(tractInfo.getBBox()).getCorners()
    skyCoords = tractWcs.pixelToSky(tractCorners)
    patchInfoDict["tract"] = (skyCoords, np.nan)

    return patchInfoDict


def generateSummaryStatsVisit(cat, colName, visitSummaryTable):
    """Generate a summary statistic in each patch or detector.

    Parameters
    ----------
    cat : `pandas.core.frame.DataFrame`
        A dataframe of the data to be plotted.
    colName : `str`
        The name of the column to be plotted.
    visitSummaryTable : `pandas.core.frame.DataFrame`
        A dataframe of the visit summary table.

    Returns
    -------
    visitInfoDict : `dict`
        A dictionary of the visit information.
    """
    visitInfoDict = {}
    for ccd in cat.detector.unique():
        if ccd is None:
            continue
        onCcd = cat["detector"] == ccd
        stat = nanMedian(cat[colName].values[onCcd])

        sumRow = visitSummaryTable["id"] == ccd
        corners = zip(visitSummaryTable["raCorners"][sumRow][0], visitSummaryTable["decCorners"][sumRow][0])
        cornersOut = []
        for ra, dec in corners:
            corner = SpherePoint(ra, dec, units=degrees)
            cornersOut.append(corner)

        visitInfoDict[ccd] = (cornersOut, stat)

    return visitInfoDict


# Inspired by matplotlib.testing.remove_ticks_and_titles
def get_and_remove_axis_text(ax) -> Tuple[List[str], List[np.ndarray]]:
    """Remove text from an Axis and its children and return with line points.

    Parameters
    ----------
    ax : `plt.Axis`
        A matplotlib figure axis.

    Returns
    -------
    texts : `List[str]`
        A list of all text strings (title and axis/legend/tick labels).
    line_xys : `List[numpy.ndarray]`
        A list of all line ``_xy`` attributes (arrays of shape ``(N, 2)``).
    """
    line_xys = [line._xy for line in ax.lines]
    texts = [text.get_text() for text in (ax.title, ax.xaxis.label, ax.yaxis.label)]
    ax.set_title("")
    ax.set_xlabel("")
    ax.set_ylabel("")

    try:
        texts_legend = ax.get_legend().texts
        texts.extend(text.get_text() for text in texts_legend)
        for text in texts_legend:
            text.set_alpha(0)
    except AttributeError:
        pass

    for idx in range(len(ax.texts)):
        texts.append(ax.texts[idx].get_text())
        ax.texts[idx].set_text("")

    ax.xaxis.set_major_formatter(null_formatter)
    ax.xaxis.set_minor_formatter(null_formatter)
    ax.yaxis.set_major_formatter(null_formatter)
    ax.yaxis.set_minor_formatter(null_formatter)
    try:
        ax.zaxis.set_major_formatter(null_formatter)
        ax.zaxis.set_minor_formatter(null_formatter)
    except AttributeError:
        pass
    for child in ax.child_axes:
        texts_child, lines_child = get_and_remove_axis_text(child)
        texts.extend(texts_child)

    return texts, line_xys


def get_and_remove_figure_text(figure: Figure):
    """Remove text from a Figure and its Axes and return with line points.

    Parameters
    ----------
    figure : `matplotlib.pyplot.Figure`
        A matplotlib figure.

    Returns
    -------
    texts : `List[str]`
        A list of all text strings (title and axis/legend/tick labels).
    line_xys : `List[numpy.ndarray]`, (N, 2)
        A list of all line ``_xy`` attributes (arrays of shape ``(N, 2)``).
    """
    texts = [str(figure._suptitle)]
    lines = []
    figure.suptitle("")

    texts.extend(text.get_text() for text in figure.texts)
    figure.texts = []

    for ax in figure.get_axes():
        texts_ax, lines_ax = get_and_remove_axis_text(ax)
        texts.extend(texts_ax)
        lines.extend(lines_ax)

    return texts, lines


def parsePlotInfo(plotInfo: Mapping[str, str]) -> str:
    """Extract information from the plotInfo dictionary and parses it into
    a meaningful string that can be added to a figure.

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
    photocalibDataset = "None"
    astroDataset = "None"

    run = plotInfo["run"]
    datasetsUsed = f"\nPhotoCalib: {photocalibDataset}, Astrometry: {astroDataset}"
    tableType = f"\nTable: {plotInfo['tableName']}"

    dataIdText = ""
    if "tract" in plotInfo.keys():
        dataIdText += f", Tract: {plotInfo['tract']}"
    if "visit" in plotInfo.keys():
        dataIdText += f", Visit: {plotInfo['visit']}"

    bandText = ""
    for band in plotInfo["bands"]:
        bandText += band + ", "
    bandsText = f", Bands: {bandText[:-2]}"
    infoText = f"\n{run}{datasetsUsed}{tableType}{dataIdText}{bandsText}"

    # Find S/N and mag keys, if present.
    snKeys = []
    magKeys = []
    selectionKeys = []
    selectionPrefix = "Selection: "
    for key, value in plotInfo.items():
        if "SN" in key or "S/N" in key:
            snKeys.append(key)
        elif "Mag" in key:
            magKeys.append(key)
        elif key.startswith(selectionPrefix):
            selectionKeys.append(key)
    # Add S/N and mag values to label, if present.
    # TODO: Do something if there are multiple sn/mag keys. Log? Warn?
    newline = "\n"
    if snKeys:
        infoText = f"{infoText}{newline if magKeys else ', '}{snKeys[0]}{plotInfo.get(snKeys[0])}"
    if magKeys:
        infoText = f"{infoText}, {magKeys[0]}{plotInfo.get(magKeys[0])}"
    if selectionKeys:
        nPrefix = len(selectionPrefix)
        selections = ", ".join(f"{key[nPrefix:]}: {plotInfo[key]}" for key in selectionKeys)
        infoText = f"{infoText}, Selections: {selections}"

    return infoText


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


def mkColormap(colorNames):
    """Make a colormap from the list of color names.

    Parameters
    ----------
    colorNames : `list`
        A list of strings that correspond to matplotlib named colors.

    Returns
    -------
    cmap : `matplotlib.colors.LinearSegmentedColormap`
        A colormap stepping through the supplied list of names.
    """
    nums = np.linspace(0, 1, len(colorNames))
    blues = []
    greens = []
    reds = []
    for num, color in zip(nums, colorNames):
        r, g, b = colors.colorConverter.to_rgb(color)
        blues.append((num, b, b))
        greens.append((num, g, g))
        reds.append((num, r, r))

    colorDict = {"blue": blues, "red": reds, "green": greens}
    cmap = colors.LinearSegmentedColormap("newCmap", colorDict)
    return cmap


def extremaSort(xs):
    """Return the IDs of the points reordered so that those furthest from the
    median, in absolute terms, are last.

    Parameters
    ----------
    xs : `np.array`
        An array of the values to sort

    Returns
    -------
    ids : `np.array`
    """
    med = nanMedian(xs)
    dists = np.abs(xs - med)
    ids = np.argsort(dists)
    return ids


def sortAllArrays(arrsToSort, sortArrayIndex=0):
    """Sort one array and then return all the others in the associated order.

    Parameters
    ----------
    arrsToSort : `list` [`np.array`]
        A list of arrays to be simultaneously sorted based on the array in
        the list position given by ``sortArrayIndex`` (defaults to be the
        first array in the list).
    sortArrayIndex : `int`, optional
        Zero-based index indicating the array on which to base the sorting.

    Returns
    -------
    arrsToSort : `list` [`np.array`]
        The list of arrays sorted on array in list index ``sortArrayIndex``.
    """
    ids = extremaSort(arrsToSort[sortArrayIndex])
    for i, arr in enumerate(arrsToSort):
        arrsToSort[i] = arr[ids]
    return arrsToSort


def addSummaryPlot(fig, loc, sumStats, label):
    """Add a summary subplot to the figure.

    Parameters
    ----------
    fig : `matplotlib.figure.Figure`
        The figure that the summary plot is to be added to.
    loc : `matplotlib.gridspec.SubplotSpec` or `int` or `(int, int, index`
        Describes the location in the figure to put the summary plot,
        can be a gridspec SubplotSpec, a 3 digit integer where the first
        digit is the number of rows, the second is the number of columns
        and the third is the index. This is the same for the tuple
        of int, int, index.
    sumStats : `dict`
        A dictionary where the patchIds are the keys which store the R.A.
        and the dec of the corners of the patch, along with a summary
        statistic for each patch.
    label : `str`
        The label to be used for the colorbar.

    Returns
    -------
    fig : `matplotlib.figure.Figure`
    """
    # Add the subplot to the relevant place in the figure
    # and sort the axis out
    axCorner = fig.add_subplot(loc)
    axCorner.yaxis.tick_right()
    axCorner.yaxis.set_label_position("right")
    axCorner.xaxis.tick_top()
    axCorner.xaxis.set_label_position("top")
    axCorner.set_aspect("equal")

    # Plot the corners of the patches and make the color
    # coded rectangles for each patch, the colors show
    # the median of the given value in the patch
    patches = []
    colors = []
    for dataId in sumStats.keys():
        (corners, stat) = sumStats[dataId]
        ra = corners[0][0].asDegrees()
        dec = corners[0][1].asDegrees()
        xy = (ra, dec)
        width = corners[2][0].asDegrees() - ra
        height = corners[2][1].asDegrees() - dec
        patches.append(Rectangle(xy, width, height))
        colors.append(stat)
        ras = [ra.asDegrees() for (ra, dec) in corners]
        decs = [dec.asDegrees() for (ra, dec) in corners]
        axCorner.plot(ras + [ras[0]], decs + [decs[0]], "k", lw=0.5)
        cenX = ra + width / 2
        cenY = dec + height / 2
        if dataId != "tract":
            axCorner.annotate(dataId, (cenX, cenY), color="k", fontsize=4, ha="center", va="center")

    # Set the bad color to transparent and make a masked array
    cmapPatch = cm.coolwarm.copy()
    cmapPatch.set_bad(color="none")
    colors = np.ma.array(colors, mask=np.isnan(colors))
    collection = PatchCollection(patches, cmap=cmapPatch)
    collection.set_array(colors)
    axCorner.add_collection(collection)

    # Add some labels
    axCorner.set_xlabel("R.A. (deg)", fontsize=7)
    axCorner.set_ylabel("Dec. (deg)", fontsize=7)
    axCorner.tick_params(axis="both", labelsize=6, length=0, pad=1.5)
    axCorner.invert_xaxis()

    # Add a colorbar
    pos = axCorner.get_position()
    yOffset = (pos.y1 - pos.y0) / 3
    cax = fig.add_axes([pos.x0, pos.y1 + yOffset, pos.x1 - pos.x0, 0.025])
    fig.colorbar(collection, cax=cax, orientation="horizontal")
    cax.text(
        0.5,
        0.48,
        label,
        color="k",
        transform=cax.transAxes,
        rotation="horizontal",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=6,
    )
    cax.tick_params(
        axis="x", labelsize=6, labeltop=True, labelbottom=False, bottom=False, top=True, pad=0.5, length=2
    )

    return fig


def shorten_list(numbers: Iterable[int], *, range_indicator: str = "-", range_separator: str = ",") -> str:
    """Shorten an iterable of integers.

    Parameters
    ----------
    numbers : `~collections.abc.Iterable` [`int`]
        Any iterable (list, set, tuple, numpy.array) of integers.
    range_indicator : `str`, optional
        The string to use to indicate a range of numbers.
    range_separator : `str`, optional
        The string to use to separate ranges of numbers.

    Returns
    -------
    result : `str`
        A shortened string representation of the list.

    Examples
    --------
    >>> shorten_list([1,2,3,5,6,8])
    "1-3,5-6,8"

    >>> shorten_list((1,2,3,5,6,8,9,10,11), range_separator=", ")
    "1-3, 5-6, 8-11"

    >>> shorten_list(range(4), range_indicator="..")
    "0..3"
    """
    # Sort the list in ascending order.
    numbers = sorted(numbers)

    if not numbers:  # empty container
        return ""

    # Initialize an empty list to hold the results to be returned.
    result = []

    # Initialize variables to track the current start and end of a list.
    start = 0
    end = 0  # initialize to 0 to handle single element lists.

    # Iterate through the sorted list of numbers
    for end in range(1, len(numbers)):
        # If the current number is the same or consecutive to the previous
        # number, skip to the next iteration.
        if numbers[end] > numbers[end - 1] + 1:  # > is used to handle duplicates, if any.
            # If the current number is not consecutive to the previous number,
            # add the current range to the result and reset the start to end.
            if start == end - 1:
                result.append(str(numbers[start]))
            else:
                result.append(range_indicator.join((str(numbers[start]), str(numbers[end - 1]))))

            # Update start.
            start = end

    # Add the final range to the result.
    if start == end:
        result.append(str(numbers[start]))
    else:
        result.append(range_indicator.join((str(numbers[start]), str(numbers[end]))))

    # Return the shortened string representation.
    return range_separator.join(result)


class PanelConfig(Config):
    """Configuration options for the plot panels used by DiaSkyPlot.

    The defaults will produce a good-looking single panel plot.
    The subplot2grid* fields correspond to matplotlib.pyplot.subplot2grid.
    """

    topSpinesVisible = Field[bool](
        doc="Draw line and ticks on top of panel?",
        default=False,
    )
    bottomSpinesVisible = Field[bool](
        doc="Draw line and ticks on bottom of panel?",
        default=True,
    )
    leftSpinesVisible = Field[bool](
        doc="Draw line and ticks on left side of panel?",
        default=True,
    )
    rightSpinesVisible = Field[bool](
        doc="Draw line and ticks on right side of panel?",
        default=True,
    )
    subplot2gridShapeRow = Field[int](
        doc="Number of rows of the grid in which to place axis.",
        default=10,
    )
    subplot2gridShapeColumn = Field[int](
        doc="Number of columns of the grid in which to place axis.",
        default=10,
    )
    subplot2gridLocRow = Field[int](
        doc="Row of the axis location within the grid.",
        default=1,
    )
    subplot2gridLocColumn = Field[int](
        doc="Column of the axis location within the grid.",
        default=1,
    )
    subplot2gridRowspan = Field[int](
        doc="Number of rows for the axis to span downwards.",
        default=5,
    )
    subplot2gridColspan = Field[int](
        doc="Number of rows for the axis to span to the right.",
        default=5,
    )


def plotProjectionWithBinning(
    ax,
    xs,
    ys,
    zs,
    cmap,
    xMin,
    xMax,
    yMin,
    yMax,
    xNumBins=45,
    yNumBins=None,
    fixAroundZero=False,
    nPointBinThresh=5000,
    isSorted=False,
    vmin=None,
    vmax=None,
    showExtremeOutliers=True,
    scatPtSize=7,
):
    """Plot color-mapped data in projection and with binning when appropriate.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes`
        Axis on which to plot the projection data.
    xs, ys : `np.array`
        Arrays containing the x and y positions of the data.
    zs : `np.array`
        Array containing the scaling value associated with the (``xs``, ``ys``)
        positions.
    cmap : `matplotlib.colors.Colormap`
        Colormap for the ``zs`` values.
    xMin, xMax, yMin, yMax : `float`
        Data limits within which to compute bin sizes.
    xNumBins : `int`, optional
        The number of bins along the x-axis.
    yNumBins : `int`, optional
        The number of bins along the y-axis. If `None`, this is set to equal
        ``xNumBins``.
    nPointBinThresh : `int`, optional
        Threshold number of points above which binning will be implemented
        for the plotting. If the number of data points is lower than this
        threshold, a basic scatter plot will be generated.
    isSorted : `bool`, optional
        Whether the data have been sorted in ``zs`` (the sorting is to
        accommodate the overplotting of points in the upper and lower
        extrema of the data).
    vmin, vmax : `float`, optional
        The min and max limits for the colorbar.
    showExtremeOutliers: `bool`, default True
        Use overlaid scatter points to show the x-y positions of the 15%
        most extreme values.
    scatPtSize : `float`, optional
        The point size to use if just plotting a regular scatter plot.

    Returns
    -------
    plotOut : `matplotlib.collections.PathCollection`
        The plot object with ``ax`` updated with data plotted here.
    """
    med = nanMedian(zs)
    mad = nanSigmaMad(zs)
    if vmin is None:
        vmin = med - 2 * mad
    if vmax is None:
        vmax = med + 2 * mad
    if fixAroundZero:
        scaleEnd = np.max([np.abs(vmin), np.abs(vmax)])
        vmin = -1 * scaleEnd
        vmax = scaleEnd

    yNumBins = xNumBins if yNumBins is None else yNumBins

    xBinEdges = np.linspace(xMin, xMax, xNumBins + 1)
    yBinEdges = np.linspace(yMin, yMax, yNumBins + 1)
    finiteMask = np.isfinite(zs)
    xs = xs[finiteMask]
    ys = ys[finiteMask]
    zs = zs[finiteMask]
    binnedStats, xEdges, yEdges, binNums = binned_statistic_2d(
        xs, ys, zs, statistic="median", bins=(xBinEdges, yBinEdges)
    )

    if len(xs) >= nPointBinThresh:
        s = min(10, max(0.5, nPointBinThresh / 10 / (len(xs) ** 0.5)))
        lw = (s**0.5) / 10
        plotOut = ax.imshow(
            binnedStats.T,
            cmap=cmap,
            extent=[xEdges[0], xEdges[-1], yEdges[-1], yEdges[0]],
            vmin=vmin,
            vmax=vmax,
        )
        if not isSorted:
            sortedArrays = sortAllArrays([zs, xs, ys])
            zs, xs, ys = sortedArrays[0], sortedArrays[1], sortedArrays[2]
        if len(xs) > 1:
            if showExtremeOutliers:
                # Find the most extreme 15% of points. The list is ordered
                # by the distance from the median, this is just the
                # head/tail 15% of points.
                extremes = int(np.floor((len(xs) / 100)) * 85)
                plotOut = ax.scatter(
                    xs[extremes:],
                    ys[extremes:],
                    c=zs[extremes:],
                    s=s,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    edgecolor="white",
                    linewidths=lw,
                )
    else:
        plotOut = ax.scatter(
            xs,
            ys,
            c=zs,
            cmap=cmap,
            s=scatPtSize,
            vmin=vmin,
            vmax=vmax,
            edgecolor="white",
            linewidths=0.2,
        )
    return plotOut
