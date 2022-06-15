import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import median_absolute_deviation as sigmaMad
import pandas as pd
from sklearn.neighbors import KernelDensity
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pathEffects

import lsst.pipe.base as pipeBase
from lsst.pex.config.pexConfig import Field, ListField, DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.pipe.tasks.dataFrameActions import MagColumnNanoJansky

from .calcFunctors import ExtinctionCorrectedMagDiff
from . import dataSelectors as dataSelectors
from .plotUtils import parsePlotInfo, addPlotInfo, stellarLocusFit, perpDistance, mkColormap


class ColorColorFitPlot(PlotAction):

    xAxisLabel = Field(doc="Label to use for the x axis", dtype=str, optional=False)
    yAxisLabel = Field(doc="Label to use for the y axis", dtype=str, optional=False)
    magLabel = Field(doc="Label to use for the magnitudes used to color code by", dtype=str, optional=False)

    plotTypes = ListField(
        doc="Selection of types of objects to plot. Can take any combination of"
            " stars, galaxies, unknown, mag, any.",
        dtype=str,
        optional=False,
        itemCheck=_validatePlotTypes,
    )

    stellarLocusFitDict = DictField(
        doc="The parameters to use for the stellar locus fit. The default parameters are examples and are "
            "not useful for any of the fits. The dict needs to contain xMin/xMax/yMin/yMax which are the "
            "limits of the initial box for fitting the stellar locus, mHW and bHW are the initial "
            "intercept and gradient for the fitting.",
        keytype=str,
        itemtype=float,
        default={"xMin": 0.1, "xMax": 0.2, "yMin": 0.1, "yMax": 0.2, "mHW": 0.5, "bHW": 0.0}
    )


    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        base.append(("x", Vector))
        base.append(("y", Vector))
        base.append(("mag", Vector))
        base.append(("snThreshold", Scalar))
        base.append(("snFlux", Vector))
        base.append((""))

        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """ NOTE currently can only check that something is not a scalar, not
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

    def makePlot(self, data: keyedData, plotInfo: Optional[Mapping[str, str]] = None,
                 fitParams: Optional[Mapping[str, float]] = None, **kwargs) -> Figure:
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
        axHist = fig. add_axes([0.65, 0.51, 0.3, 0.31])
        xs = data["x"]
        ys = data["y"]
        mags = data["mag"]

        # TODO: Make a no data fig function and use here
        if len(xs) == 0 or len(ys) == 0:
            return fig

        # Points to use for the fit
        fitPoints = np.where((xs > fitParams["xMin"]) & (xs < fitParams["xMax"])
                             & (ys > fitParams["yMin"]) & (ys < fitParams["yMax"]))[0]

        # Plot the initial fit box
        ax.plot([fitParams["xMin"], fitParams["xMax"], fitParams["xMax"], fitParams["xMin"],
                fitParams["xMin"]], [fitParams["yMin"], fitParams["yMin"], fitParams["yMax"],
                fitParams["yMax"], fitParams["yMin"]], "k", alpha=0.3)

        # Add some useful information to the plot
        bbox = dict(alpha=0.9, facecolor="white", edgecolor="none")
        medMag = np.median(mags)

        try:
            ids = (data["fluxSn"] < data["snThreshold"]*1.1)
            medMag = np.nanmedian(mags[ids])
        except AttributeError:
            SN = "NA"
            medMag = "NA"

        infoText = "N Used: {}\nN Total: {}\nS/N cut: {}\n".format(len(fitPoints), len(data["x"]), SN)
        infoText += r"Mag $\lesssim$: " + "{:0.2f}".format(medMag)
        ax.text(0.05, 0.78, infoText, color="k", transform=ax.transAxes,
                fontsize=8, bbox=bbox)

        # Calculate the density of the points
        xy = np.vstack([xs, ys]).T
        kde = KernelDensity(kernel="gaussian").fit(xy)
        z = np.exp(kde.score_samples(xy))

        ax.scatter(xs[~fitPoints], ys[~fitPoints], c=z[~fitPoints], cmap="binary", s=0.3)
        fitScatter = ax.scatter(xs[fitPoints], ys[fitPoints], c=z[fitPoints], cmap=newBlues,
                                label="Used for Fit", s=0.3)

        # Add colorbar
        cbAx = fig.add_axes([0.12, 0.08, 0.43, 0.04])
        plt.colorbar(fitScatter, cax=cbAx, orientation="horizontal")
        cbText = cbAx.text(0.5, 0.5, "Number Density", color="k", rotation="horizontal",
                           transform=cbAx.transAxes, ha="center", va="center", fontsize=8)
        cbText.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])
        cbAx.set_xticks([np.min(z[fitPoints]), np.max(z[fitPoints])], labels=["Less", "More"])

        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)

        # Set useful axis limits
        percsX = np.nanpercentile(xs, [0.5, 99.5])
        percsY = np.nanpercentile(ys, [0.5, 99.5])
        x5 = (percsX[1] - percsX[0])/5
        y5 = (percsY[1] - percsY[0])/5
        ax.set_xlim(percsX[0] - x5, percsX[1] + x5)
        ax.set_ylim(percsY[0] - y5, percsY[1] + y5)

        # Plot the fit lines
        if np.fabs(fitParams["mHW"]) > 1:
            ysFitLineHW = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLineHW = (ysFitLineHW - fitParams["bHW"])/fitParams["mHW"]
            ysFitLine = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLine = (ysFitLine - fitParams["bODR"])/fitParams["mODR"]
            ysFitLine2 = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLine2 = (ysFitLine2 - fitParams["bODR2"])/fitParams["mODR2"]

        else:
            xsFitLineHW = np.array([fitParams["xMin"], fitParams["xMax"]])
            ysFitLineHW = fitParams["mHW"]*xsFitLineHW + fitParams["bHW"]
            xsFitLine = [fitParams["xMin"], fitParams["xMax"]]
            ysFitLine = [fitParams["mODR"]*xsFitLine[0] + fitParams["bODR"],
                         fitParams["mODR"]*xsFitLine[1] + fitParams["bODR"]]
            xsFitLine2 = [fitParams["xMin"], fitParams["xMax"]]
            ysFitLine2 = [fitParams["mODR2"]*xsFitLine2[0] + fitParams["bODR2"],
                          fitParams["mODR2"]*xsFitLine2[1] + fitParams["bODR2"]]

        ax.plot(xsFitLineHW, ysFitLineHW, "w", lw=2)
        lineHW, = ax.plot(xsFitLineHW, ysFitLineHW, "g", lw=1, ls="--", label="Hardwired")

        ax.plot(xsFitLine, ysFitLine, "w", lw=2)
        lineInit, = ax.plot(xsFitLine, ysFitLine, "b", lw=1, ls="--", label="Initial")

        ax.plot(xsFitLine2, ysFitLine2, "w", lw=2)
        lineRefit, = ax.plot(xsFitLine2, ysFitLine2, "k", lw=1, ls="--", label="Refit")

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
        if np.fabs(fitParams["mHW"]) > 1:
            xMid = (fitParams["yMin"] - fitParams["bODR2"])/fitParams["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"]*xs + fitParams["bPerpMin"]
        else:
            xs = np.array([fitParams["xMin"] - 0.2, fitParams["xMin"], fitParams["xMin"] + 0.2])
            ys = xs*fitParams["mPerp"] + fitParams["bPerpMin"]
        ax.plot(xs, ys, "k--", alpha=0.7)

        if np.fabs(fitParams["mHW"]) > 1:
            xMid = (fitParams["yMax"] - fitParams["bODR2"])/fitParams["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"]*xs + fitParams["bPerpMax"]
        else:
            xs = np.array([fitParams["xMax"] - 0.2, fitParams["xMax"], fitParams["xMax"] + 0.2])
            ys = xs*fitParams["mPerp"] + fitParams["bPerpMax"]
        ax.plot(xs, ys, "k--", alpha=0.7)

        # Add a histogram
        axHist.set_ylabel("Number")
        axHist.set_xlabel("Distance to Line Fit")
        medDists = np.median(dists)
        madDists = sigmaMad(dists)
        meanDists = np.mean(dists)

        rmsDists = np.sqrt(np.mean(np.array(dists)**2))
        axHist.set_xlim(meanDists - 2.0*rmsDists, meanDists + 2.0*rmsDists)
        lineMedian = axHist.axvline(medDists, color="k", label="Median: {:0.3f}".format(medDists))
        lineMad = axHist.axvline(medDists + madDists, color="k", ls="--",
                                 label="sigma MAD: {:0.3f}".format(madDists))
        axHist.axvline(medDists - madDists, color="k", ls="--")
        lineMean = axHist.axvline(meanDists, color="C0", label="Mean: {:0.3f}".format(meanDists))
        lineRms = axHist.axvline(meanDists + rmsDists, color="C0", ls="--",
                                 label="RMS: {:0.3f}".format(rmsDists))
        axHist.axvline(meanDists - rmsDists, color="C0", ls="--")

        linesForLegend = [lineHW, lineInit, lineRefit, fitScatter, lineMedian, lineMad, lineMean, lineRms]
        fig.legend(handles=linesForLegend, fontsize=8, bbox_to_anchor=(1.0, 0.99),
                   bbox_transform=fig.transFigure, ncol=2)

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
        minXs = -1*np.min(np.fabs(percsDists))
        maxXs = np.min(np.fabs(percsDists))
        plotPoints = ((dists < maxXs) & (dists > minXs))
        xs = np.array(dists)[plotPoints]
        ys = mags[fitPoints][plotPoints]
        H, xEdges, yEdges = np.histogram2d(xs, ys, bins=(11, 11))
        xBinWidth = xEdges[1] - xEdges[0]
        yBinWidth = yEdges[1] - yEdges[0]
        axContour.contour(xEdges[:-1] + xBinWidth/2, yEdges[:-1] + yBinWidth/2, H.T, levels=7,
                          cmap=newBlues)
        axContour.set_xlabel("Distance to Line Fit")
        axContour.set_ylabel(self.magLabel)

        fig = addPlotInfo(plt.gcf(), plotInfo)

        return fig
