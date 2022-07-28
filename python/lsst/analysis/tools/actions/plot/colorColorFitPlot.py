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

from typing import Mapping, Optional, cast

import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
from lsst.analysis.tools import PlotAction
from lsst.pex.config import Field, ListField
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from scipy.stats import median_absolute_deviation as sigmaMad
from sklearn.neighbors import KernelDensity

from ...interfaces import KeyedData, KeyedDataSchema, Scalar, Vector
from .plotUtils import addPlotInfo, mkColormap, perpDistance


class ColorColorFitPlot(PlotAction):

    xAxisLabel = Field[str](doc="Label to use for the x axis", optional=False)
    yAxisLabel = Field[str](doc="Label to use for the y axis", optional=False)
    magLabel = Field[str](doc="Label to use for the magnitudes used to color code by", optional=False)

    plotTypes = ListField[str](
        doc="Selection of types of objects to plot. Can take any combination of"
        " stars, galaxies, unknown, mag, any.",
        optional=False,
    )

    plotName = Field[str](doc="The name for the plot.", optional=False)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | type[Scalar]]] = []
        base.append(("x", Vector))
        base.append(("y", Vector))
        base.append(("mag", Vector))
        base.append(("approxMagDepth", Scalar))
        base.append((f"{self.plotName}_sigmaMAD", Scalar))
        base.append((f"{self.plotName}_median", Scalar))
        base.append((f"{self.plotName}_hardwired_sigmaMAD", Scalar))
        base.append((f"{self.plotName}_hardwired_median", Scalar))
        base.append(("xMin", Scalar))
        base.append(("xMax", Scalar))
        base.append(("yMin", Scalar))
        base.append(("yMax", Scalar))
        base.append(("mHW", Scalar))
        base.append(("bHW", Scalar))
        base.append(("mODR", Scalar))
        base.append(("bODR", Scalar))
        base.append(("yBoxMin", Scalar))
        base.append(("yBoxMax", Scalar))
        base.append(("bPerpMin", Scalar))
        base.append(("bPerpMax", Scalar))
        base.append(("mODR2", Scalar))
        base.append(("bODR2", Scalar))
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
        plotInfo: Optional[Mapping[str, str]] = None,
        **kwargs,
    ) -> Figure:
        """Make stellar locus plots using pre fitted values.

        Parameters
        ----------
        catPlot : `pandas.core.frame.DataFrame`
            The catalog to plot the points from.
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:
                ``"run"``
                    The output run for the plots (`str`).
                ``"skymap"``
                    The type of skymap used for the data (`str`).
                ``"filter"``
                    The filter used for this data (`str`).
                ``"tract"``
                    The tract that the data comes from (`str`).
        fitParams : `dict`
            The parameters of the fit to the stellar locus calculated
            elsewhere, they are used to plot the fit line on the
            figure.
                ``"bHW"``
                    The hardwired intercept to fall back on.
                ``"b_odr"``
                    The intercept calculated by the orthogonal distance
                    regression fitting.
                ``"mHW"``
                    The hardwired gradient to fall back on.
                ``"m_odr"``
                    The gradient calculated by the orthogonal distance
                    regression fitting.
                ``"magLim"``
                    The magnitude limit used in the fitting.
                ``"x1`"``
                    The x minimum of the box used in the fit.
                ``"x2"``
                    The x maximum of the box used in the fit.
                ``"y1"``
                    The y minimum of the box used in the fit.
                ``"y2"``
                    The y maximum of the box used in the fit.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.

        Notes
        -----
        Makes a color-color plot of `self.config.xColName` against
        `self.config.yColName`, these points are color coded by i band
        CModel magnitude. The stellar locus fits calculated from
        the calcStellarLocus task are then overplotted. The axis labels
        are given by `self.config.xLabel` and `self.config.yLabel`.
        The selector given in `self.config.sourceSelectorActions`
        is used for source selection. The distance of the points to
        the fit line is given in a histogram in the second panel.
        """

        # Define a new colormap
        newBlues = mkColormap(["paleturquoise", "midnightblue"])

        # Make a figure with three panels
        fig = plt.figure(dpi=300)
        ax = fig.add_axes([0.12, 0.25, 0.43, 0.60])
        axContour = fig.add_axes([0.65, 0.11, 0.3, 0.31])
        axHist = fig.add_axes([0.65, 0.51, 0.3, 0.31])
        xs = cast(Vector, data["x"])
        ys = cast(Vector, data["y"])
        mags = data["mag"]

        # TODO: Make a no data fig function and use here
        if len(xs) == 0 or len(ys) == 0:
            return fig

        # Points to use for the fit
        # type ignore because Vector needs a prototype interface
        fitPoints = np.where(
            (xs > data["xMin"])  # type: ignore
            & (xs < data["xMax"])  # type: ignore
            & (ys > data["yMin"])  # type: ignore
            & (ys < data["yMax"])  # type: ignore
        )[0]

        # Plot the initial fit box
        ax.plot(
            [data["xMin"], data["xMax"], data["xMax"], data["xMin"], data["xMin"]],
            [data["yMin"], data["yMin"], data["yMax"], data["yMax"], data["yMin"]],
            "k",
            alpha=0.3,
        )

        # Add some useful information to the plot
        bbox = dict(alpha=0.9, facecolor="white", edgecolor="none")
        medMag = np.median(cast(Vector, mags))

        # TODO: GET THE SN FROM THE EARLIER PREP STEP
        SN = "-"
        infoText = "N Used: {}\nN Total: {}\nS/N cut: {}\n".format(len(fitPoints), len(xs), SN)
        infoText += r"Mag $\lesssim$: " + "{:0.2f}".format(medMag)
        ax.text(0.05, 0.78, infoText, color="k", transform=ax.transAxes, fontsize=8, bbox=bbox)

        # Calculate the density of the points
        xy = np.vstack([xs, ys]).T
        kde = KernelDensity(kernel="gaussian").fit(xy)
        z = np.exp(kde.score_samples(xy))

        ax.scatter(xs[~fitPoints], ys[~fitPoints], c=z[~fitPoints], cmap="binary", s=0.3)
        fitScatter = ax.scatter(
            xs[fitPoints], ys[fitPoints], c=z[fitPoints], cmap=newBlues, label="Used for Fit", s=0.3
        )

        # Add colorbar
        cbAx = fig.add_axes([0.12, 0.08, 0.43, 0.04])
        plt.colorbar(fitScatter, cax=cbAx, orientation="horizontal")
        cbText = cbAx.text(
            0.5,
            0.5,
            "Number Density",
            color="k",
            rotation="horizontal",
            transform=cbAx.transAxes,
            ha="center",
            va="center",
            fontsize=8,
        )
        cbText.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])
        cbAx.set_xticks([np.min(z[fitPoints]), np.max(z[fitPoints])], labels=["Less", "More"])

        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)

        # Set useful axis limits
        percsX = np.nanpercentile(xs, [0.5, 99.5])
        percsY = np.nanpercentile(ys, [0.5, 99.5])
        x5 = (percsX[1] - percsX[0]) / 5
        y5 = (percsY[1] - percsY[0]) / 5
        ax.set_xlim(percsX[0] - x5, percsX[1] + x5)
        ax.set_ylim(percsY[0] - y5, percsY[1] + y5)

        # Plot the fit lines
        if np.fabs(data["mHW"]) > 1:
            ysFitLineHW = np.array([data["yMin"], data["yMax"]])
            xsFitLineHW = (ysFitLineHW - data["bHW"]) / data["mHW"]
            ysFitLine = np.array([data["yMin"], data["yMax"]])
            xsFitLine = (ysFitLine - data["bODR"]) / data["mODR"]
            ysFitLine2 = np.array([data["yMin"], data["yMax"]])
            xsFitLine2 = (ysFitLine2 - data["bODR2"]) / data["mODR2"]

        else:
            xsFitLineHW = np.array([data["xMin"], data["xMax"]])
            ysFitLineHW = data["mHW"] * xsFitLineHW + data["bHW"]  # type: ignore
            xsFitLine = np.array([data["xMin"], data["xMax"]])
            ysFitLine = np.array(
                [
                    data["mODR"] * xsFitLine[0] + data["bODR"],
                    data["mODR"] * xsFitLine[1] + data["bODR"],
                ]
            )
            xsFitLine2 = np.array([data["xMin"], data["xMax"]])
            ysFitLine2 = np.array(
                [
                    data["mODR2"] * xsFitLine2[0] + data["bODR2"],
                    data["mODR2"] * xsFitLine2[1] + data["bODR2"],
                ]
            )

        ax.plot(xsFitLineHW, ysFitLineHW, "w", lw=2)
        (lineHW,) = ax.plot(xsFitLineHW, ysFitLineHW, "g", lw=1, ls="--", label="Hardwired")

        ax.plot(xsFitLine, ysFitLine, "w", lw=2)
        (lineInit,) = ax.plot(xsFitLine, ysFitLine, "b", lw=1, ls="--", label="Initial")

        ax.plot(xsFitLine2, ysFitLine2, "w", lw=2)
        (lineRefit,) = ax.plot(xsFitLine2, ysFitLine2, "k", lw=1, ls="--", label="Refit")

        # Calculate the distances to that line
        # Need two points to characterise the lines we want
        # to get the distances to
        p1 = np.array([xsFitLine[0], ysFitLine[0]])
        p2 = np.array([xsFitLine[1], ysFitLine[1]])

        p1HW = np.array([xsFitLine[0], ysFitLineHW[0]])
        p2HW = np.array([xsFitLine[1], ysFitLineHW[1]])

        distsHW = perpDistance(p1HW, p2HW, zip(xs[fitPoints], ys[fitPoints]))
        dists = perpDistance(p1, p2, zip(xs[fitPoints], ys[fitPoints]))

        # Now we have the information for the perpendicular line we
        # can use it to calculate the points at the ends of the
        # perpendicular lines that intersect at the box edges
        if np.fabs(data["mHW"]) > 1:
            xMid = (data["yMin"] - data["bODR2"]) / data["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = data["mPerp"] * xs + data["bPerpMin"]
        else:
            xs = np.array([data["xMin"] - 0.2, data["xMin"], data["xMin"] + 0.2])
            ys = xs * data["mPerp"] + data["bPerpMin"]
        ax.plot(xs, ys, "k--", alpha=0.7)

        if np.fabs(data["mHW"]) > 1:
            xMid = (data["yMax"] - data["bODR2"]) / data["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = data["mPerp"] * xs + data["bPerpMax"]
        else:
            xs = np.array([data["xMax"] - 0.2, data["xMax"], data["xMax"] + 0.2])
            ys = xs * data["mPerp"] + data["bPerpMax"]
        ax.plot(xs, ys, "k--", alpha=0.7)

        # Add a histogram
        axHist.set_ylabel("Number")
        axHist.set_xlabel("Distance to Line Fit")
        medDists = np.median(dists)
        madDists = sigmaMad(dists)
        meanDists = np.mean(dists)

        axHist.set_xlim(meanDists - 2.0 * madDists, meanDists + 2.0 * madDists)
        lineMedian = axHist.axvline(medDists, color="k", label="Median: {:0.3f}".format(medDists))
        lineMad = axHist.axvline(
            medDists + madDists, color="k", ls="--", label="sigma MAD: {:0.3f}".format(madDists)
        )
        axHist.axvline(medDists - madDists, color="k", ls="--")

        linesForLegend = [lineHW, lineInit, lineRefit, fitScatter, lineMedian, lineMad]
        fig.legend(
            handles=linesForLegend,
            fontsize=8,
            bbox_to_anchor=(1.0, 0.99),
            bbox_transform=fig.transFigure,
            ncol=2,
        )

        axHist.hist(dists, bins=100, histtype="step", label="Refit", color="C0")
        axHist.hist(distsHW, bins=100, histtype="step", label="HW", color="C0", alpha=0.5)

        alphas = [1.0, 0.5]
        handles = [Rectangle((0, 0), 1, 1, color="none", ec="C0", alpha=a) for a in alphas]
        labels = ["Refit", "HW"]
        axHist.legend(handles, labels, fontsize=6, loc="upper right")

        # Add a contour plot showing the magnitude dependance
        # of the distance to the fit
        axContour.invert_yaxis()
        axContour.axvline(0.0, color="k", ls="--", zorder=-1)
        percsDists = np.nanpercentile(dists, [4, 96])
        minXs = -1 * np.min(np.fabs(percsDists))
        maxXs = np.min(np.fabs(percsDists))
        plotPoints = (dists < maxXs) & (dists > minXs)
        xs = np.array(dists)[plotPoints]
        ys = cast(Vector, cast(Vector, mags)[cast(Vector, fitPoints)])[cast(Vector, plotPoints)]
        H, xEdges, yEdges = np.histogram2d(xs, ys, bins=(11, 11))
        xBinWidth = xEdges[1] - xEdges[0]
        yBinWidth = yEdges[1] - yEdges[0]
        axContour.contour(
            xEdges[:-1] + xBinWidth / 2, yEdges[:-1] + yBinWidth / 2, H.T, levels=7, cmap=newBlues
        )
        axContour.set_xlabel("Distance to Line Fit")
        axContour.set_ylabel(self.magLabel)

        fig = addPlotInfo(plt.gcf(), plotInfo)

        return fig
