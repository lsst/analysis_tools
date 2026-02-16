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

__all__ = ("ColorColorFitPlot",)

from collections.abc import Mapping
from typing import cast

import matplotlib.patheffects as pathEffects
import numpy as np
import scipy.stats
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from scipy.ndimage import median_filter

from lsst.pex.config import Field, ListField, RangeField
from lsst.utils.plotting import make_figure, set_rubin_plotstyle

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from ...math import nanMean, nanMedian, nanSigmaMad
from ..keyedData.stellarLocusFit import perpDistance
from .plotUtils import addPlotInfo, mkColormap


class ColorColorFitPlot(PlotAction):
    """Make a color-color plot and overplot a prefited line to the fit region.

    This is mostly used for the stellar locus plots and also includes panels
    that illustrate the goodness of the given fit.
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis", optional=False)
    yAxisLabel = Field[str](doc="Label to use for the y axis", optional=False)
    magLabel = Field[str](doc="Label to use for the magnitudes used to color code by", optional=False)

    plotTypes = ListField[str](
        doc="Selection of types of objects to plot. Can take any combination of"
        " stars, galaxies, unknown, mag, any.",
        default=["stars"],
    )

    plotName = Field[str](doc="The name for the plot.", optional=False)
    minPointsForFit = RangeField[int](
        doc="Minimum number of valid objects to bother attempting a fit.",
        default=5,
        min=1,
        deprecated="This field is no longer used. The value should go as an "
        "entry to the paramsDict keyed as minObjectForFit.  Will be removed "
        "after v27.",
    )

    xLims = ListField[float](
        doc="Minimum and maximum x-axis limit to force (provided as a list of [xMin, xMax]). "
        "If `None`, limits will be computed and set based on the data.",
        dtype=float,
        default=None,
        optional=True,
    )

    yLims = ListField[float](
        doc="Minimum and maximum y-axis limit to force (provided as a list of [yMin, yMax]). "
        "If `None`, limits will be computed and set based on the data.",
        dtype=float,
        default=None,
        optional=True,
    )

    doPlotRedBlueHists = Field[bool](
        doc="Plot distance from fit histograms separated into blue and red star subsamples?",
        default=False,
        optional=True,
    )

    doPlotDistVsColor = Field[bool](
        doc="Plot distance from fit as a function of color in lower right panel?",
        default=True,
        optional=True,
    )

    publicationStyle = Field[bool](
        doc="Use a publication-quality plot style.",
        default=False,
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | type[Scalar]]] = []
        base.append(("x", Vector))
        base.append(("y", Vector))
        base.append(("mag", Vector))
        base.append(("approxMagDepth", Scalar))
        base.append((f"{self.plotName}_sigmaMAD", Scalar))
        base.append((f"{self.plotName}_median", Scalar))
        base.append(("mODR", Scalar))
        base.append(("bODR", Scalar))
        base.append(("bPerpMin", Scalar))
        base.append(("bPerpMax", Scalar))
        base.append(("mPerp", Scalar))

        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """NOTE currently can only check that something is not a scalar, not
        check that data is consistent with Vector
        """
        needed = self.getInputSchema(**kwargs)
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not in input")
        for name, typ in needed:
            isScalar = issubclass((colType := type(data[name.format(**kwargs)])), Scalar)
            if isScalar and typ != Scalar:
                raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}")

    def makePlot(
        self,
        data: KeyedData,
        plotInfo: Mapping[str, str],
        **kwargs,
    ) -> Figure:
        """Make stellar locus plots using pre fitted values.

        Parameters
        ----------
        data : `KeyedData`
            The data to plot the points from, for more information
            please see the notes section.
        plotInfo : `dict`
            A dictionary of information about the data being plotted
            with keys:

            * ``"run"``
                  The output run for the plots (`str`).
            * ``"skymap"``
                  The type of skymap used for the data (`str`).
            * ``"filter"``
                  The filter used for this data (`str`).
            * ``"tract"``
                  The tract that the data comes from (`str`).

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.

        Notes
        -----
        The axis labels are given by `self.xAxisLabel` and `self.yAxisLabel`.
        The perpendicular distance of the points to the fit line is given in a
        histogram in the second panel.

        For the code to work it expects various quantities to be present in
        the `data` that it is given.

        The quantities that are expected to be present are:

         * Statistics that are shown on the plot or used by the plotting code:
            * ``approxMagDepth``
                  The approximate magnitude corresponding to the SN cut used.
            * ``f"{self.plotName}_sigmaMAD"``
                  The sigma mad of the distances to the line fit.
            * ``f"{self.plotName or ''}_median"``
                  The median of the distances to the line fit.

         * Parameters from the fitting code that are illustrated on the plot:
            * ``"bFixed"``
                  The fixed intercept to fall back on.
            * ``"mFixed"``
                  The fixed gradient to fall back on.
            * ``"bODR"``
                  The intercept calculated by the final orthogonal distance
                  regression fitting.
            * ``"mODR"``
                  The gradient calculated by the final orthogonal distance
                  regression fitting.
            * ``"xMin`"``
                  The x minimum of the box used in the fit.
            * ``"xMax"``
                  The x maximum of the box used in the fit.
            * ``"yMin"``
                  The y minimum of the box used in the fit.
            * ``"yMax"``
                  The y maximum of the box used in the fit.
            * ``"mPerp"``
                  The gradient of the line perpendicular to the line from
                  the second ODR fit.
            * ``"bPerpMin"``
                  The intercept of the perpendicular line that goes through
                  xMin.
            * ``"bPerpMax"``
                  The intercept of the perpendicular line that goes through
                  xMax.
            * ``"goodPoints"``
                  The points that passed the initial set of cuts (typically
                  in fluxType S/N, extendedness, magnitude, and isfinite).
            * ``"fitPoints"``
                  The points use in the final fit.

          * The main inputs to plot:
                x, y, mag

        Examples
        --------
        An example of the plot produced from this code is here:

        .. image:: /_static/analysis_tools/stellarLocusExample.png

        For a detailed example of how to make a plot from the command line
        please see the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """
        set_rubin_plotstyle()
        paramDict = data.pop("paramDict")
        # Points to use for the fit.
        fitPoints = data.pop("fitPoints")
        # Points with finite values for x, y, and mag.
        goodPoints = data.pop("goodPoints")

        # TODO: Make a no data fig function and use here.
        if sum(fitPoints) < paramDict["minObjectForFit"]:
            fig = make_figure(dpi=120)
            ax = fig.add_axes([0.12, 0.25, 0.43, 0.62])
            ax.tick_params(labelsize=7)
            noDataText = (
                "Number of objects after cuts ({})\nis less than the minimum required\nby "
                "paramDict[minObjectForFit] ({})".format(sum(fitPoints), int(paramDict["minObjectForFit"]))
            )
            fig.text(0.5, 0.5, noDataText, ha="center", va="center", fontsize=8)
            fig = addPlotInfo(fig, plotInfo)
            return fig

        # Define new colormaps.
        newBlues = mkColormap(["#3D5195"])
        newGrays = mkColormap(["#CFCFCF", "#595959"])

        # Make a figure with three panels.
        fig = make_figure(dpi=300)
        ax = fig.add_axes([0.12, 0.25, 0.43, 0.62])
        if self.doPlotDistVsColor:
            axLowerRight = fig.add_axes([0.65, 0.11, 0.26, 0.34])
        else:
            axLowerRight = fig.add_axes([0.65, 0.11, 0.3, 0.34])
        axHist = fig.add_axes([0.65, 0.55, 0.3, 0.32])

        xs = cast(Vector, data["x"])
        ys = cast(Vector, data["y"])
        mags = cast(Vector, data["mag"])

        # Plot the initial fit box.
        (initialBox,) = ax.plot(
            [paramDict["xMin"], paramDict["xMax"], paramDict["xMax"], paramDict["xMin"], paramDict["xMin"]],
            [paramDict["yMin"], paramDict["yMin"], paramDict["yMax"], paramDict["yMax"], paramDict["yMin"]],
            "k",
            alpha=0.3,
            label="Initial selection",
        )

        if not self.publicationStyle:
            # Add some useful information to the plot.
            bbox = dict(alpha=0.9, facecolor="white", edgecolor="none")
            infoText = f"N Total: {sum(goodPoints)}\nN Used: {sum(fitPoints)}"
            ax.text(0.04, 0.97, infoText, color="k", transform=ax.transAxes, fontsize=7, bbox=bbox, va="top")

        # Calculate the point density for the Used and NotUsed subsamples.
        xyUsed = np.vstack([xs[fitPoints], ys[fitPoints]])
        xyNotUsed = np.vstack([xs[~fitPoints & goodPoints], ys[~fitPoints & goodPoints]])

        # Try using a Gaussian KDE to get color bars showing number density. If
        # there are not enough points for KDE, use a constant color.
        zUsed = np.ones(np.sum(fitPoints))
        if len(xyUsed.ravel()) >= 2:
            try:
                zUsed = scipy.stats.gaussian_kde(xyUsed)(xyUsed)
            except np.linalg.LinAlgError:
                pass

        zNotUsed = np.ones(np.sum(~fitPoints & goodPoints))
        if len(xyNotUsed.ravel()) >= 2:
            try:
                zNotUsed = scipy.stats.gaussian_kde(xyNotUsed)(xyNotUsed)
            except np.linalg.LinAlgError:
                pass

        notUsedScatter = ax.scatter(
            xs[~fitPoints & goodPoints], ys[~fitPoints & goodPoints], c=zNotUsed, cmap=newGrays, s=0.3
        )
        fitScatter = ax.scatter(
            xs[fitPoints], ys[fitPoints], c=zUsed, cmap=newBlues, s=0.3, label="Used for Fit"
        )

        # Add colorbars.
        cbAx = fig.add_axes([0.12, 0.07, 0.43, 0.04])
        fig.colorbar(fitScatter, cax=cbAx, orientation="horizontal")
        cbKwargs = {
            "color": "k",
            "rotation": "horizontal",
            "ha": "center",
            "va": "center",
            "fontsize": 7,
        }
        cbText = cbAx.text(
            0.5,
            0.5,
            "Number Density (used in fit)",
            transform=cbAx.transAxes,
            **cbKwargs,
        )
        cbText.set_path_effects([pathEffects.Stroke(linewidth=1.5, foreground="w"), pathEffects.Normal()])
        cbAx.set_xticks([np.min(zUsed), np.max(zUsed)], labels=["Less", "More"], fontsize=7)
        cbAxNotUsed = fig.add_axes([0.12, 0.11, 0.43, 0.04])
        fig.colorbar(notUsedScatter, cax=cbAxNotUsed, orientation="horizontal")
        cbText = cbAxNotUsed.text(
            0.5,
            0.5,
            "Number Density (not used in fit)",
            transform=cbAxNotUsed.transAxes,
            **cbKwargs,
        )
        cbText.set_path_effects([pathEffects.Stroke(linewidth=1.5, foreground="w"), pathEffects.Normal()])
        cbAxNotUsed.set_xticks([])

        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)
        ax.tick_params(labelsize=7)

        # Set axis limits from configs if set, otherwise based on the data.
        if self.xLims is not None:
            ax.set_xlim(self.xLims[0], self.xLims[1])
        else:
            percsX = np.nanpercentile(xs[goodPoints], [0.5, 99.5])
            x5 = (percsX[1] - percsX[0]) / 5
            ax.set_xlim(percsX[0] - x5, percsX[1] + x5)
        if self.yLims is not None:
            ax.set_ylim(self.yLims[0], self.yLims[1])
        else:
            percsY = np.nanpercentile(ys[goodPoints], [0.5, 99.5])
            y5 = (percsY[1] - percsY[0]) / 5
            ax.set_ylim(percsY[0] - y5, percsY[1] + y5)

        # Plot the fit lines.
        if np.fabs(paramDict["mFixed"]) > 1:
            ysFitLineFixed = np.array([paramDict["yMin"], paramDict["yMax"]])
            xsFitLineFixed = (ysFitLineFixed - paramDict["bFixed"]) / paramDict["mFixed"]
            ysFitLine = np.array([paramDict["yMin"], paramDict["yMax"]])
            xsFitLine = (ysFitLine - data["bODR"]) / data["mODR"]

        else:
            xsFitLineFixed = np.array([paramDict["xMin"], paramDict["xMax"]])
            ysFitLineFixed = paramDict["mFixed"] * xsFitLineFixed + paramDict["bFixed"]
            xsFitLine = np.array([paramDict["xMin"], paramDict["xMax"]])
            ysFitLine = np.array(
                [
                    data["mODR"] * xsFitLine[0] + data["bODR"],
                    data["mODR"] * xsFitLine[1] + data["bODR"],
                ]
            )

        ax.plot(xsFitLineFixed, ysFitLineFixed, "w", lw=1.5)
        (lineFixed,) = ax.plot(xsFitLineFixed, ysFitLineFixed, "#C85200", lw=1, ls="--", label="Fixed")
        ax.plot(xsFitLine, ysFitLine, "w", lw=1.5)
        (lineOdrFit,) = ax.plot(xsFitLine, ysFitLine, "k", lw=1, ls="--", label="ODR Fit")
        ax.legend(
            handles=[initialBox, lineFixed, lineOdrFit], handlelength=1.5, fontsize=6, loc="lower right"
        )

        # Calculate the distances (in mmag) to the line for the data used in
        # the fit. Two points are needed to characterize the lines we want
        # to get the distances to.
        p1 = np.array([xsFitLine[0], ysFitLine[0]])
        p2 = np.array([xsFitLine[1], ysFitLine[1]])

        p1Fixed = np.array([xsFitLineFixed[0], ysFitLineFixed[0]])
        p2Fixed = np.array([xsFitLineFixed[1], ysFitLineFixed[1]])

        # Convert to millimags.
        statsUnitStr = "mmag"
        distsFixed = np.array(perpDistance(p1Fixed, p2Fixed, zip(xs[fitPoints], ys[fitPoints]))) * 1000
        dists = np.array(perpDistance(p1, p2, zip(xs[fitPoints], ys[fitPoints]))) * 1000
        maxDist = np.abs(np.nanmax(dists)) / 1000  # These will be used to set the fit boundary line limits.
        minDist = np.abs(np.nanmin(dists)) / 1000
        # Now we have the information for the perpendicular line we can use it
        # to calculate the points at the ends of the perpendicular lines that
        # intersect at the box edges.
        if np.fabs(paramDict["mFixed"]) > 1:
            xMid = (paramDict["yMin"] - data["bODR"]) / data["mODR"]
            xsFit = np.array([xMid - max(0.2, maxDist), xMid, xMid + max(0.2, minDist)])
            ysFit = data["mPerp"] * xsFit + data["bPerpMin"]
        else:
            xsFit = np.array(
                [
                    paramDict["xMin"] - max(0.2, np.fabs(paramDict["mFixed"]) * maxDist),
                    paramDict["xMin"],
                    paramDict["xMin"] + max(0.2, np.fabs(paramDict["mFixed"]) * minDist),
                ]
            )
            ysFit = xsFit * data["mPerp"] + data["bPerpMin"]
        ax.plot(xsFit, ysFit, "k--", alpha=0.7, lw=1)

        if np.fabs(paramDict["mFixed"]) > 1:
            xMid = (paramDict["yMax"] - data["bODR"]) / data["mODR"]
            xsFit = np.array([xMid - max(0.2, maxDist), xMid, xMid + max(0.2, minDist)])
            ysFit = data["mPerp"] * xsFit + data["bPerpMax"]
        else:
            xsFit = np.array(
                [
                    paramDict["xMax"] - max(0.2, np.fabs(paramDict["mFixed"]) * maxDist),
                    paramDict["xMax"],
                    paramDict["xMax"] + max(0.2, np.fabs(paramDict["mFixed"]) * minDist),
                ]
            )
            ysFit = xsFit * data["mPerp"] + data["bPerpMax"]
        ax.plot(xsFit, ysFit, "k--", alpha=0.7, lw=1)

        # Compute statistics for fit.
        medDists = nanMedian(dists)
        madDists = nanSigmaMad(dists)
        meanDists = nanMean(dists)
        rmsDists = np.sqrt(np.mean(np.array(dists) ** 2))

        xMid = paramDict["xMin"] + 0.5 * (paramDict["xMax"] - paramDict["xMin"])
        if self.doPlotRedBlueHists:
            blueStars = (xs[fitPoints] < xMid) & (xs[fitPoints] >= paramDict["xMin"])
            blueDists = dists[blueStars]
            blueMedDists = nanMedian(blueDists)
            redStars = (xs[fitPoints] >= xMid) & (xs[fitPoints] <= paramDict["xMax"])
            redDists = dists[redStars]
            redMedDists = nanMedian(redDists)

        if self.doPlotRedBlueHists:
            blueStars = (xs[fitPoints] < xMid) & (xs[fitPoints] >= paramDict["xMin"])
            blueDists = dists[blueStars]
            blueMedDists = nanMedian(blueDists)
            redStars = (xs[fitPoints] >= xMid) & (xs[fitPoints] <= paramDict["xMax"])
            redDists = dists[redStars]
            redMedDists = nanMedian(redDists)

        # Add a histogram.
        axHist.set_ylabel("Number", fontsize=7)
        axHist.set_xlabel(f"Distance to Line Fit ({statsUnitStr})", fontsize=7)
        axHist.tick_params(labelsize=7)
        nSigToPlot = 3.5
        axHist.set_xlim(meanDists - nSigToPlot * madDists, meanDists + nSigToPlot * madDists)
        lineMedian = axHist.axvline(
            medDists, color="k", lw=1, alpha=0.5, label=f"Median: {medDists:0.2f} {statsUnitStr}"
        )
        lineMad = axHist.axvline(
            medDists + madDists,
            color="k",
            ls="--",
            lw=1,
            alpha=0.5,
            label=r"$\sigma_{MAD}$" + f": {madDists:0.2f} {statsUnitStr}",
        )
        axHist.axvline(medDists - madDists, color="k", ls="--", lw=1, alpha=0.5)
        lineRms = axHist.axvline(
            meanDists + rmsDists,
            color="k",
            ls=":",
            lw=1,
            alpha=0.3,
            label=f"RMS: {rmsDists:0.2f} {statsUnitStr}",
        )
        axHist.axvline(meanDists - rmsDists, color="k", ls=":", lw=1, alpha=0.3)
        if self.doPlotRedBlueHists:
            lineBlueMedian = axHist.axvline(
                blueMedDists,
                color="blue",
                ls=":",
                lw=1.0,
                alpha=0.5,
                label=f"blueMed: {blueMedDists:0.2f}",
            )
            lineRedMedian = axHist.axvline(
                redMedDists,
                color="red",
                ls=":",
                lw=1.0,
                alpha=0.5,
                label=f"redMed: {redMedDists:0.2f}",
            )
            linesForLegend = [lineMedian, lineMad, lineRms, lineBlueMedian, lineRedMedian]
        else:
            linesForLegend = [lineMedian, lineMad, lineRms]
        fig.legend(
            handles=linesForLegend,
            handlelength=1.0,
            fontsize=6,
            loc="lower right",
            bbox_to_anchor=(0.955, 0.89),
            bbox_transform=fig.transFigure,
            ncol=2,
        )

        axHist.hist(dists, bins=100, histtype="stepfilled", label="ODR Fit", color="k", ec="k", alpha=0.3)
        axHist.hist(distsFixed, bins=100, histtype="step", label="Fixed", color="#C85200", alpha=1.0)
        if self.doPlotRedBlueHists:
            axHist.hist(blueDists, bins=100, histtype="stepfilled", color="blue", ec="blue", alpha=0.3)
            axHist.hist(redDists, bins=100, histtype="stepfilled", color="red", ec="red", alpha=0.3)

        handles = [Rectangle((0, 0), 1, 1, color="k", alpha=0.4)]
        handles.append(Rectangle((0, 0), 1, 1, color="none", ec="#C85200", alpha=1.0))
        labels = ["ODR Fit", "Fixed"]
        if self.doPlotRedBlueHists:
            handles.append(Rectangle((0, 0), 1, 1, color="blue", alpha=0.3))
            handles.append(Rectangle((0, 0), 1, 1, color="red", alpha=0.3))
            labels = ["ODR Fit", "Blue Stars", "Red Stars", "Fixed"]
        axHist.legend(handles, labels, fontsize=5, loc="upper right")

        if self.doPlotDistVsColor:
            axLowerRight.axhline(0.0, color="k", ls="-", lw=0.8, zorder=-1)
            # Compute and plot a running median of dists vs. color.
            if np.fabs(data["mODR"]) > 1.0:
                xRun = ys[fitPoints].copy()
                axLowerRight.set_xlabel(self.yAxisLabel, fontsize=7)
            else:
                xRun = xs[fitPoints].copy()
                axLowerRight.set_xlabel(self.xAxisLabel, fontsize=7)
            lowerRightPlot = axLowerRight.scatter(xRun, dists, c=mags[fitPoints], cmap=newBlues, s=0.2)
            yRun = dists.copy()
            xySorted = zip(xRun, yRun)
            xySorted = sorted(xySorted)
            xSorted = [x for x, y in xySorted]
            ySorted = [y for x, y in xySorted]
            nCumulate = int(max(3, len(xRun) // 10))
            yRunMedian = median_filter(ySorted, size=nCumulate)
            axLowerRight.plot(xSorted, yRunMedian, "w", lw=1.8)
            axLowerRight.plot(xSorted, yRunMedian, c="#595959", ls="-", lw=1.1, label="Running Median")
            axLowerRight.set_ylim(-2.5 * madDists, 2.5 * madDists)
            axLowerRight.set_ylabel(f"Distance to Line Fit ({statsUnitStr})", fontsize=7)
            axLowerRight.legend(fontsize=4, loc="upper right", handlelength=1.0)
            # Add colorbars.
            cbAx = fig.add_axes([0.915, 0.11, 0.014, 0.34])
            fig.colorbar(lowerRightPlot, cax=cbAx, orientation="vertical")
            cbKwargs = {
                "color": "k",
                "rotation": "vertical",
                "ha": "center",
                "va": "center",
                "fontsize": 4,
            }
            cbText = cbAx.text(
                0.5,
                0.5,
                self.magLabel,
                transform=cbAx.transAxes,
                **cbKwargs,
            )
            cbText.set_path_effects([pathEffects.Stroke(linewidth=1.2, foreground="w"), pathEffects.Normal()])
            cbAx.tick_params(length=2, labelsize=3.5)
        else:
            # Add a contour plot showing the magnitude dependance of the
            # distance to the fit.
            axLowerRight.invert_yaxis()
            axLowerRight.axvline(0.0, color="k", ls="--", zorder=-1)
            percsDists = np.nanpercentile(dists, [4, 96])
            minXs = -1 * np.min(np.fabs(percsDists))
            maxXs = np.min(np.fabs(percsDists))
            plotPoints = (dists < maxXs) & (dists > minXs)
            xsContour = np.array(dists)[plotPoints]
            ysContour = cast(Vector, cast(Vector, mags)[cast(Vector, fitPoints)])[cast(Vector, plotPoints)]
            H, xEdges, yEdges = np.histogram2d(xsContour, ysContour, bins=(11, 11))
            xBinWidth = xEdges[1] - xEdges[0]
            yBinWidth = yEdges[1] - yEdges[0]
            axLowerRight.contour(
                xEdges[:-1] + xBinWidth / 2, yEdges[:-1] + yBinWidth / 2, H.T, levels=7, cmap=newBlues
            )
            axLowerRight.set_xlabel(f"Distance to Line Fit ({statsUnitStr})", fontsize=8)
            axLowerRight.set_ylabel(self.magLabel, fontsize=8)
            axLowerRight.set_xlim(meanDists - nSigToPlot * madDists, meanDists + nSigToPlot * madDists)
        axLowerRight.tick_params(labelsize=6)

        # This is here because matplotlib occasionally decides
        # that there needs to be 10^15 minor ticks.
        # This is probably unneeded and may make the plot look
        # rather busy if it ever finds enough memory to render.
        # If this ever becomes a problem it can be looked at in
        # the future.
        for ax in fig.get_axes():
            ax.minorticks_off()

        fig.canvas.draw()
        if not self.publicationStyle:
            fig = addPlotInfo(fig, plotInfo)

        return fig
