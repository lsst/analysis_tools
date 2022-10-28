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

from typing import List, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.odr as scipyODR
from lsst.geom import Box2D, SpherePoint, degrees
from lsst.pex.config import Config, Field
from matplotlib import colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

null_formatter = matplotlib.ticker.NullFormatter()


def parsePlotInfo(dataId, runName, tableName, bands, plotName, SN):
    """Parse plot info from the dataId
    Parameters
    ----------
    dataId : `lsst.daf.butler.core.dimensions.`
             `_coordinate._ExpandedTupleDataCoordinate`
    runName : `str`
    Returns
    -------
    plotInfo : `dict`
    """
    plotInfo = {"run": runName, "tableName": tableName, "plotName": plotName, "SN": SN}

    for dataInfo in dataId:
        plotInfo[dataInfo.name] = dataId[dataInfo.name]

    bandStr = ""
    for band in bands:
        bandStr += ", " + band
    plotInfo["bands"] = bandStr[2:]

    if "tract" not in plotInfo.keys():
        plotInfo["tract"] = "N/A"
    if "visit" not in plotInfo.keys():
        plotInfo["visit"] = "N/A"

    return plotInfo


def generateSummaryStats(data, skymap, plotInfo):
    """Generate a summary statistic in each patch or detector
    Parameters
    ----------
    data : `dict`
    skymap : `lsst.skymap.ringsSkyMap.RingsSkyMap`
    plotInfo : `dict`

    Returns
    -------
    patchInfoDict : `dict`
    """

    # TODO: what is the more generic type of skymap?
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
        stat = np.nanmedian(data[yCol][onPatch])
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


def generateSummaryStatsVisit(cat, colName, visitSummaryTable, plotInfo):
    """Generate a summary statistic in each patch or detector
    Parameters
    ----------
    cat : `pandas.core.frame.DataFrame`
    colName : `str`
    visitSummaryTable : `pandas.core.frame.DataFrame`
    plotInfo : `dict`
    Returns
    -------
    visitInfoDict : `dict`
    """

    visitInfoDict = {}
    for ccd in cat.detector.unique():
        if ccd is None:
            continue
        onCcd = cat["detector"] == ccd
        stat = np.nanmedian(cat[colName].values[onCcd])

        sumRow = visitSummaryTable["id"] == ccd
        corners = zip(visitSummaryTable["raCorners"][sumRow][0], visitSummaryTable["decCorners"][sumRow][0])
        cornersOut = []
        for (ra, dec) in corners:
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


def get_and_remove_figure_text(figure: plt.Figure):
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


def addPlotInfo(fig, plotInfo):
    """Add useful information to the plot
    Parameters
    ----------
    fig : `matplotlib.figure.Figure`
    plotInfo : `dict`
    Returns
    -------
    fig : `matplotlib.figure.Figure`
    """

    # TO DO: figure out how to get this information
    photocalibDataset = "None"
    astroDataset = "None"

    fig.text(0.01, 0.99, plotInfo["plotName"], fontsize=8, transform=fig.transFigure, ha="left", va="top")

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
    SNText = f", S/N: {plotInfo.get('SN', 'N/A')}"
    infoText = f"\n{run}{datasetsUsed}{tableType}{dataIdText}{bandsText}{SNText}"
    fig.text(0.01, 0.98, infoText, fontsize=7, transform=fig.transFigure, alpha=0.6, ha="left", va="top")

    return fig


def stellarLocusFit(xs, ys, paramDict):
    """Make a fit to the stellar locus
    Parameters
    ----------
    xs : `numpy.ndarray`
        The color on the xaxis
    ys : `numpy.ndarray`
        The color on the yaxis
    paramDict : lsst.pex.config.dictField.Dict
        A dictionary of parameters for line fitting
        xMin : `float`
            The minimum x edge of the box to use for initial fitting
        xMax : `float`
            The maximum x edge of the box to use for initial fitting
        yMin : `float`
            The minimum y edge of the box to use for initial fitting
        yMax : `float`
            The maximum y edge of the box to use for initial fitting
        mHW : `float`
            The hardwired gradient for the fit
        bHW : `float`
            The hardwired intercept of the fit
    Returns
    -------
    paramsOut : `dict`
        A dictionary of the calculated fit parameters
        xMin : `float`
            The minimum x edge of the box to use for initial fitting
        xMax : `float`
            The maximum x edge of the box to use for initial fitting
        yMin : `float`
            The minimum y edge of the box to use for initial fitting
        yMax : `float`
            The maximum y edge of the box to use for initial fitting
        mHW : `float`
            The hardwired gradient for the fit
        bHW : `float`
            The hardwired intercept of the fit
        mODR : `float`
            The gradient calculated by the ODR fit
        bODR : `float`
            The intercept calculated by the ODR fit
        yBoxMin : `float`
            The y value of the fitted line at xMin
        yBoxMax : `float`
            The y value of the fitted line at xMax
        bPerpMin : `float`
            The intercept of the perpendicular line that goes through xMin
        bPerpMax : `float`
            The intercept of the perpendicular line that goes through xMax
        mODR2 : `float`
            The gradient from the second round of fitting
        bODR2 : `float`
            The intercept from the second round of fitting
        mPerp : `float`
            The gradient of the line perpendicular to the line from the
            second fit
    Notes
    -----
    The code does two rounds of fitting, the first is initiated using the
    hardwired values given in the `paramDict` parameter and is done using
    an Orthogonal Distance Regression fit to the points defined by the
    box of xMin, xMax, yMin and yMax. Once this fitting has been done a
    perpendicular bisector is calculated at either end of the line and
    only points that fall within these lines are used to recalculate the fit.
    """

    # Points to use for the fit
    fitPoints = np.where(
        (xs > paramDict["xMin"])
        & (xs < paramDict["xMax"])
        & (ys > paramDict["yMin"])
        & (ys < paramDict["yMax"])
    )[0]

    linear = scipyODR.polynomial(1)

    data = scipyODR.Data(xs[fitPoints], ys[fitPoints])
    odr = scipyODR.ODR(data, linear, beta0=[paramDict["bHW"], paramDict["mHW"]])
    params = odr.run()
    mODR = float(params.beta[1])
    bODR = float(params.beta[0])

    paramsOut = {
        "xMin": paramDict["xMin"],
        "xMax": paramDict["xMax"],
        "yMin": paramDict["yMin"],
        "yMax": paramDict["yMax"],
        "mHW": paramDict["mHW"],
        "bHW": paramDict["bHW"],
        "mODR": mODR,
        "bODR": bODR,
    }

    # Having found the initial fit calculate perpendicular ends
    mPerp = -1.0 / mODR
    # When the gradient is really steep we need to use
    # the y limits of the box rather than the x ones

    if np.abs(mODR) > 1:
        yBoxMin = paramDict["yMin"]
        xBoxMin = (yBoxMin - bODR) / mODR
        yBoxMax = paramDict["yMax"]
        xBoxMax = (yBoxMax - bODR) / mODR
    else:
        yBoxMin = mODR * paramDict["xMin"] + bODR
        xBoxMin = paramDict["xMin"]
        yBoxMax = mODR * paramDict["xMax"] + bODR
        xBoxMax = paramDict["xMax"]

    bPerpMin = yBoxMin - mPerp * xBoxMin

    paramsOut["yBoxMin"] = yBoxMin
    paramsOut["bPerpMin"] = bPerpMin

    bPerpMax = yBoxMax - mPerp * xBoxMax

    paramsOut["yBoxMax"] = yBoxMax
    paramsOut["bPerpMax"] = bPerpMax

    # Use these perpendicular lines to chose the data and refit
    fitPoints = (ys > mPerp * xs + bPerpMin) & (ys < mPerp * xs + bPerpMax)
    data = scipyODR.Data(xs[fitPoints], ys[fitPoints])
    odr = scipyODR.ODR(data, linear, beta0=[bODR, mODR])
    params = odr.run()
    mODR = float(params.beta[1])
    bODR = float(params.beta[0])

    paramsOut["mODR2"] = float(params.beta[1])
    paramsOut["bODR2"] = float(params.beta[0])

    paramsOut["mPerp"] = -1.0 / paramsOut["mODR2"]

    return paramsOut


def perpDistance(p1, p2, points):
    """Calculate the perpendicular distance to a line from a point
    Parameters
    ----------
    p1 : `numpy.ndarray`
        A point on the line
    p2 : `numpy.ndarray`
        Another point on the line
    points : `zip`
        The points to calculate the distance to
    Returns
    -------
    dists : `list`
        The distances from the line to the points. Uses the cross
        product to work this out.
    """
    dists = []
    for point in points:
        point = np.array(point)
        distToLine = np.cross(p2 - p1, point - p1) / np.linalg.norm(p2 - p1)
        dists.append(distToLine)

    return dists


def mkColormap(colorNames):
    """Make a colormap from the list of color names.
    Parameters
    ----------
    colorNames : `list`
        A list of strings that correspond to matplotlib
        named colors.
    Returns
    -------
    cmap : `matplotlib.colors.LinearSegmentedColormap`
    """

    nums = np.linspace(0, 1, len(colorNames))
    blues = []
    greens = []
    reds = []
    for (num, color) in zip(nums, colorNames):
        r, g, b = colors.colorConverter.to_rgb(color)
        blues.append((num, b, b))
        greens.append((num, g, g))
        reds.append((num, r, r))

    colorDict = {"blue": blues, "red": reds, "green": greens}
    cmap = colors.LinearSegmentedColormap("newCmap", colorDict)
    return cmap


def extremaSort(xs):
    """Return the ids of the points reordered so that those
    furthest from the median, in absolute terms, are last.
    Parameters
    ----------
    xs : `np.array`
        An array of the values to sort
    Returns
    -------
    ids : `np.array`
    """

    med = np.median(xs)
    dists = np.abs(xs - med)
    ids = np.argsort(dists)
    return ids


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
    cmapPatch = plt.cm.coolwarm.copy()
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
    plt.colorbar(collection, cax=cax, orientation="horizontal")
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
