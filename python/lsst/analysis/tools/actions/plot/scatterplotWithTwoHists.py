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

from functools import partial
from itertools import chain
from typing import Mapping, NamedTuple, Optional, cast

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sps
from lsst.analysis.tools.actions.scalar.scalarActions import CountAction, MedianAction, SigmaMadAction
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField
from matplotlib import gridspec
from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from matplotlib.figure import Figure
from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ...interfaces import (
    KeyedData,
    KeyedDataAction,
    KeyedDataSchema,
    PlotAction,
    Scalar,
    ScalarAction,
    Vector,
    VectorAction,
)
from ..keyedData import KeyedScalars
from ..vector import SnSelector
from .plotUtils import addPlotInfo, mkColormap

# from .plotUtils import addSummaryPlot, generateSummaryStats

# ignore because coolwarm is actually part of module
cmapPatch = plt.cm.coolwarm.copy()  # type: ignore
cmapPatch.set_bad(color="none")

sigmaMad = partial(sps.median_abs_deviation, scale="normal")  # type: ignore


class _ApproxMedian(ScalarAction):
    vectorKey = Field[str](doc="Key for the vector to perform action on", optional=False)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.vectorKey.format(**kwargs), Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])
        x = int(len(value) / 10)
        return np.nanmedian(value[-x:])


class _StatsImpl(KeyedScalars):
    vectorKey = Field[str](doc="Column key to compute scalars")

    snFluxType = Field[str](doc="column key for the flux type used in SN selection")

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(vectorKey=self.vectorKey)
        self.scalarActions.sigmaMad = SigmaMadAction(vectorKey=self.vectorKey)
        self.scalarActions.count = CountAction(vectorKey=self.vectorKey)
        self.scalarActions.approxMag = _ApproxMedian(vectorKey=self.snFluxType)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        return super().__call__(data, **(kwargs | dict(mask=mask)))


class ScatterPlotStatsAction(KeyedDataAction):
    vectorKey = Field[str](doc="Vector on which to compute statistics")
    highSNSelector = ConfigurableActionField[VectorAction](
        doc="Selector used to determine high SN Objects", default=SnSelector(threshold=2700)
    )
    lowSNSelector = ConfigurableActionField[VectorAction](
        doc="Selector used to determine low SN Objects", default=SnSelector(threshold=500)
    )
    fluxType = Field[str](doc="Vector key to use to compute signal to noise ratio", default="{band}_psfFlux")

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield (self.vectorKey, Vector)
        yield (self.fluxType, Vector)
        yield from self.highSNSelector.getInputSchema()
        yield from self.lowSNSelector.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return (
            (f'{self.identity or ""}HighSNMask', Vector),
            (f'{self.identity or ""}LowSNMask', Vector),
            (f"{{band}}_lowSN{self.identity.capitalize() if self.identity else ''}_median", Scalar),
            (f"{{band}}_lowSN{self.identity.capitalize() if self.identity else ''}_sigmaMad", Scalar),
            (f"{{band}}_lowSN{self.identity.capitalize() if self.identity else ''}_count", Scalar),
            (f"{{band}}_lowSN{self.identity.capitalize() if self.identity else ''}_approxMag", Scalar),
            (f"{{band}}_highSN{self.identity.capitalize() if self.identity else ''}_median", Scalar),
            (f"{{band}}_highSN{self.identity.capitalize() if self.identity else ''}_sigmaMad", Scalar),
            (f"{{band}}_highSN{self.identity.capitalize() if self.identity else ''}_count", Scalar),
            (f"{{band}}_highSN{self.identity.capitalize() if self.identity else ''}_approxMag", Scalar),
            ("highThreshold", Scalar),
            ("lowThreshold", Scalar),
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        highMaskKey = f'{self.identity or ""}HighSNMask'
        results[highMaskKey] = self.highSNSelector(data, **kwargs)

        lowMaskKey = f'{self.identity or ""}LowSNMask'
        results[lowMaskKey] = self.lowSNSelector(data, **kwargs)

        prefix = f"{band}_" if (band := kwargs.get("band")) else ""

        stats = _StatsImpl(vectorKey=self.vectorKey, snFluxType=self.fluxType)
        # this is sad, but pex_config seems to have broken behavior that
        # is dangerous to fix
        stats.setDefaults()
        for maskKey, typ in ((lowMaskKey, "low"), (highMaskKey, "high")):
            for name, value in stats(data, **(kwargs | {"mask": results[maskKey]})).items():
                tmpKey = (
                    f"{prefix}{typ}SN{self.identity.capitalize() if self.identity else '' }_{name}".format(
                        **kwargs
                    )
                )
                results[tmpKey] = value
        results["highSnThreshold"] = self.highSNSelector.threshold  # type: ignore
        results["lowSnThreshold"] = self.lowSNSelector.threshold  # type: ignore

        return results


def _validatePlotTypes(value):
    return value in ("stars", "galaxies", "unknown", "any", "mag")


# ignore type because of conflicting name on tuple baseclass
class _StatsContainer(NamedTuple):
    median: Scalar
    sigmaMad: Scalar
    count: Scalar  # type: ignore
    approxMag: Scalar


class ScatterPlotWithTwoHists(PlotAction):
    yLims = ListField[float](
        doc="ylimits of the plot, if not specified determined from data",
        length=2,
        optional=True,
    )

    xLims = ListField[float](
        doc="xlimits of the plot, if not specified determined from data", length=2, optional=True
    )
    xAxisLabel = Field[str](doc="Label to use for the x axis", optional=False)
    yAxisLabel = Field[str](doc="Label to use for the y axis", optional=False)
    magLabel = Field[str](doc="Label to use for the magnitudes used for SNR", optional=False)
    nBins = Field[float](doc="Number of bins on x axis", default=40.0)
    plot2DHist = Field[bool](
        doc="Plot a 2D histogram in dense areas of points on the scatter plot."
        "Doesn't look great if plotting multiple datasets on top of each other.",
        default=True,
    )
    plotTypes = ListField[str](
        doc="Selection of types of objects to plot. Can take any combination of"
        " stars, galaxies, unknown, mag, any",
        optional=False,
        itemCheck=_validatePlotTypes,
    )

    _stats = ("median", "sigmaMad", "count", "approxMag")

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | type[Scalar]]] = []
        if "stars" in self.plotTypes:  # type: ignore
            base.append(("xStars", Vector))
            base.append(("yStars", Vector))
            base.append(("starsHighSNMask", Vector))
            base.append(("starsLowSNMask", Vector))
            # statistics
            for name in self._stats:
                base.append((f"{{band}}_highSNStars_{name}", Scalar))
                base.append((f"{{band}}_lowSNStars_{name}", Scalar))
        if "galaxies" in self.plotTypes:  # type: ignore
            base.append(("xGalaxies", Vector))
            base.append(("yGalaxies", Vector))
            base.append(("galaxiesHighSNMask", Vector))
            base.append(("galaxiesLowSNMask", Vector))
            # statistics
            for name in self._stats:
                base.append((f"{{band}}_highSNGalaxies_{name}", Scalar))
                base.append((f"{{band}}_lowSNGalaxies_{name}", Scalar))
        if "unknown" in self.plotTypes:  # type: ignore
            base.append(("xUnknown", Scalar))
            base.append(("yUnknown", Scalar))
            base.append(("unknownHighSNMask", Vector))
            base.append(("unknownLowSNMask", Vector))
            # statistics
            for name in self._stats:
                base.append((f"{{band}}_highSNUnknown_{name}", Scalar))
                base.append((f"{{band}}_lowSNUnknown_{name}", Scalar))
        if "any" in self.plotTypes:  # type: ignore
            base.extend((("x", Vector), ("y", Vector)))
            base.append(("anyHighSNMask", Vector))
            base.append(("anySNMask", Vector))
            # statistics
            for name in self._stats:
                base.append((f"{{band}}_highSNAny_{name}", Scalar))
                base.append((f"{{band}}_lowSNAny_{name}", Scalar))
        base.append(("lowSnThreshold", Scalar))
        base.append(("highSnThreshold", Scalar))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)
        # table is a dict that needs: x, y, run, skymap, filter, tract,

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """NOTE currently can only check that something is not a Scalar, not
        check that the data is consistent with Vector
        """
        needed = self.getFormattedInputSchema(**kwargs)
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in input")
        for name, typ in needed:
            isScalar = issubclass((colType := type(data[name.format(**kwargs)])), Scalar)
            if isScalar and typ != Scalar:
                raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}")

    def makePlot(
        self,
        data: KeyedData,
        plotInfo: Optional[Mapping[str, str]] = None,
        sumStats: Optional[Mapping] = None,
        **kwargs,
    ) -> Figure:
        """Makes a generic plot with a 2D histogram and collapsed histograms of
        each axis.
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
        sumStats : `dict`
            A dictionary where the patchIds are the keys which store the R.A.
            and dec of the corners of the patch, along with a summary
            statistic for each patch.
        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.
        Notes
        -----
        Uses the axisLabels config options `x` and `y` and the axisAction
        config options `xAction` and `yAction` to plot a scatter
        plot of the values against each other. A histogram of the points
        collapsed onto each axis is also plotted. A summary panel showing the
        median of the y value in each patch is shown in the upper right corner
        of the resultant plot. The code uses the selectorActions to decide
        which points to plot and the statisticSelector actions to determine
        which points to use for the printed statistics.
        """
        if not self.plotTypes:
            noDataFig = Figure()
            noDataFig.text(0.3, 0.5, "No data to plot after selectors applied")
            noDataFig = addPlotInfo(noDataFig, plotInfo)
            return noDataFig

        fig = plt.figure(dpi=300)
        gs = gridspec.GridSpec(4, 4)

        # add the various plot elements
        ax, imhist = self._scatterPlot(data, fig, gs, **kwargs)
        self._makeTopHistogram(data, fig, gs, ax, **kwargs)
        self._makeSideHistogram(data, fig, gs, ax, imhist, **kwargs)
        # Needs info from run quantum
        # sumStats = generateSummaryStats(data, key)
        # fig = addSummaryPlot(fig, gs[0, -1], sumStats, label)

        plt.draw()
        plt.subplots_adjust(wspace=0.0, hspace=0.0, bottom=0.22, left=0.21)
        fig = addPlotInfo(fig, plotInfo)
        return fig

    def _scatterPlot(
        self, data: KeyedData, fig: Figure, gs: gridspec.GridSpec, **kwargs
    ) -> tuple[Axes, Optional[PolyCollection]]:
        # Main scatter plot
        ax = fig.add_subplot(gs[1:, :-1])

        newBlues = mkColormap(["paleturquoise", "midnightBlue"])
        newReds = mkColormap(["lemonchiffon", "firebrick"])

        binThresh = 5

        yBinsOut = []
        linesForLegend = []

        toPlotList = []
        histIm = None
        highStats: _StatsContainer
        lowStats: _StatsContainer
        if "stars" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self._stats:
                highArgs[name] = cast(Scalar, data[f"{{band}}_highSNStars_{name}".format(**kwargs)])
                lowArgs[name] = cast(Scalar, data[f"{{band}}_lowSNStars_{name}".format(**kwargs)])
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["xStars"],
                    data["yStars"],
                    data["starsHighSNMask"],
                    data["starsLowSNMask"],
                    "midnightblue",
                    newBlues,
                    highStats,
                    lowStats,
                )
            )
        if "galaxies" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self._stats:
                highArgs[name] = cast(Scalar, data[f"{{band}}_highSNGalaxies_{name}".format(**kwargs)])
                lowArgs[name] = cast(Scalar, data[f"{{band}}_lowSNGalaxies_{name}".format(**kwargs)])
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["xGalaxies"],
                    data["yGalaxies"],
                    data["galaxiesHighSNMask"],
                    data["galaxiesLowSNMask"],
                    "firebrick",
                    newReds,
                    highStats,
                    lowStats,
                )
            )
        if "unknown" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self._stats:
                highArgs[name] = cast(Scalar, data[f"{{band}}_highSNUnknown_{name}".format(**kwargs)])
                lowArgs[name] = cast(Scalar, data[f"{{band}}_lowSNUnknown_{name}".format(**kwargs)])
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["xUnknown"],
                    data["yUnknown"],
                    data["unknownHighSNMask"],
                    data["unknownLowSNMask"],
                    "green",
                    None,
                    highStats,
                    lowStats,
                )
            )
        if "any" in self.plotTypes:  # type: ignore
            highArgs = {}
            lowArgs = {}
            for name in self._stats:
                highArgs[name] = cast(Scalar, data[f"{{band}}_highSNUnknown_{name}".format(**kwargs)])
                lowArgs[name] = cast(Scalar, data[f"{{band}}_lowSNUnknown_{name}".format(**kwargs)])
            highStats = _StatsContainer(**highArgs)
            lowStats = _StatsContainer(**lowArgs)

            toPlotList.append(
                (
                    data["x"],
                    data["y"],
                    data["anyHighSNMask"],
                    data["anyLowSNMask"],
                    "purple",
                    None,
                    highStats,
                    lowStats,
                )
            )

        xMin = None
        for (j, (xs, ys, highSn, lowSn, color, cmap, highStats, lowStats)) in enumerate(toPlotList):
            highSn = cast(Vector, highSn)
            lowSn = cast(Vector, lowSn)
            # ensure the columns are actually array
            xs = np.array(xs)
            ys = np.array(ys)
            sigMadYs = sigmaMad(ys, nan_policy="omit")
            if len(xs) < 2:
                (medLine,) = ax.plot(
                    xs, np.nanmedian(ys), color, label=f"Median: {np.nanmedian(ys):0.3g}", lw=0.8
                )
                linesForLegend.append(medLine)
                sigMads = np.array([sigmaMad(ys, nan_policy="omit")] * len(xs))
                (sigMadLine,) = ax.plot(
                    xs,
                    np.nanmedian(ys) + 1.0 * sigMads,
                    color,
                    alpha=0.8,
                    lw=0.8,
                    label=r"$\sigma_{MAD}$: " + f"{sigMads[0]:0.3g}",
                )
                ax.plot(xs, np.nanmedian(ys) - 1.0 * sigMads, color, alpha=0.8)
                linesForLegend.append(sigMadLine)
                histIm = None
                continue

            [xs1, xs25, xs50, xs75, xs95, xs97] = np.nanpercentile(xs, [1, 25, 50, 75, 95, 97])
            xScale = (xs97 - xs1) / 20.0  # This is ~5% of the data range

            # 40 was used as the default number of bins because it looked good
            xEdges = np.arange(
                np.nanmin(xs) - xScale,
                np.nanmax(xs) + xScale,
                (np.nanmax(xs) + xScale - (np.nanmin(xs) - xScale)) / self.nBins,
            )
            medYs = np.nanmedian(ys)
            fiveSigmaHigh = medYs + 5.0 * sigMadYs
            fiveSigmaLow = medYs - 5.0 * sigMadYs
            binSize = (fiveSigmaHigh - fiveSigmaLow) / 101.0
            yEdges = np.arange(fiveSigmaLow, fiveSigmaHigh, binSize)

            counts, xBins, yBins = np.histogram2d(xs, ys, bins=(xEdges, yEdges))
            yBinsOut.append(yBins)
            countsYs = np.sum(counts, axis=1)

            ids = np.where((countsYs > binThresh))[0]
            xEdgesPlot = xEdges[ids][1:]
            xEdges = xEdges[ids]

            if len(ids) > 1:
                # Create the codes needed to turn the sigmaMad lines
                # into a path to speed up checking which points are
                # inside the area.
                codes = np.ones(len(xEdgesPlot) * 2) * Path.LINETO
                codes[0] = Path.MOVETO
                codes[-1] = Path.CLOSEPOLY

                meds = np.zeros(len(xEdgesPlot))
                threeSigMadVerts = np.zeros((len(xEdgesPlot) * 2, 2))
                sigMads = np.zeros(len(xEdgesPlot))

                for (i, xEdge) in enumerate(xEdgesPlot):
                    ids = np.where((xs < xEdge) & (xs > xEdges[i]) & (np.isfinite(ys)))[0]
                    med = np.median(ys[ids])
                    sigMad = sigmaMad(ys[ids])
                    meds[i] = med
                    sigMads[i] = sigMad
                    threeSigMadVerts[i, :] = [xEdge, med + 3 * sigMad]
                    threeSigMadVerts[-(i + 1), :] = [xEdge, med - 3 * sigMad]

                (medLine,) = ax.plot(xEdgesPlot, meds, color, label="Running Median")
                linesForLegend.append(medLine)

                # Make path to check which points lie within one sigma mad
                threeSigMadPath = Path(threeSigMadVerts, codes)

                # Add lines for the median +/- 3 * sigma MAD
                (threeSigMadLine,) = ax.plot(
                    xEdgesPlot,
                    threeSigMadVerts[: len(xEdgesPlot), 1],
                    color,
                    alpha=0.4,
                    label=r"3$\sigma_{MAD}$",
                )
                ax.plot(xEdgesPlot[::-1], threeSigMadVerts[len(xEdgesPlot) :, 1], color, alpha=0.4)

                # Add lines for the median +/- 1 * sigma MAD
                (sigMadLine,) = ax.plot(
                    xEdgesPlot, meds + 1.0 * sigMads, color, alpha=0.8, label=r"$\sigma_{MAD}$"
                )
                linesForLegend.append(sigMadLine)
                ax.plot(xEdgesPlot, meds - 1.0 * sigMads, color, alpha=0.8)

                # Add lines for the median +/- 2 * sigma MAD
                (twoSigMadLine,) = ax.plot(
                    xEdgesPlot, meds + 2.0 * sigMads, color, alpha=0.6, label=r"2$\sigma_{MAD}$"
                )
                linesForLegend.append(twoSigMadLine)
                linesForLegend.append(threeSigMadLine)
                ax.plot(xEdgesPlot, meds - 2.0 * sigMads, color, alpha=0.6)

                # Check which points are outside 3 sigma MAD of the median
                # and plot these as points.
                inside = threeSigMadPath.contains_points(np.array([xs, ys]).T)
                ax.plot(xs[~inside], ys[~inside], ".", ms=3, alpha=0.3, mfc=color, mec=color, zorder=-1)

                # Add some stats text
                xPos = 0.65 - 0.4 * j
                bbox = dict(edgecolor=color, linestyle="--", facecolor="none")
                highThresh = data["highSnThreshold"]
                statText = f"S/N > {highThresh} Stats ({self.magLabel} < {highStats.approxMag})\n"
                highStatsStr = (
                    f"Median: {highStats.median}    "
                    + r"$\sigma_{MAD}$: "
                    + f"{highStats.sigmaMad}    "
                    + r"N$_{points}$: "
                    + f"{highStats.count}"
                )
                statText += highStatsStr
                fig.text(xPos, 0.090, statText, bbox=bbox, transform=fig.transFigure, fontsize=6)

                bbox = dict(edgecolor=color, linestyle=":", facecolor="none")
                lowThresh = data["lowSnThreshold"]
                statText = f"S/N > {lowThresh} Stats ({self.magLabel} < {lowStats.approxMag})\n"
                lowStatsStr = (
                    f"Median: {lowStats.median}    "
                    + r"$\sigma_{MAD}$: "
                    + f"{lowStats.sigmaMad}    "
                    + r"N$_{points}$: "
                    + f"{lowStats.count}"
                )
                statText += lowStatsStr
                fig.text(xPos, 0.020, statText, bbox=bbox, transform=fig.transFigure, fontsize=6)

                if self.plot2DHist:
                    histIm = ax.hexbin(xs[inside], ys[inside], gridsize=75, cmap=cmap, mincnt=1, zorder=-2)

                # If there are not many sources being used for the
                # statistics then plot them individually as just
                # plotting a line makes the statistics look wrong
                # as the magnitude estimation is iffy for low
                # numbers of sources.
                if np.sum(highSn) < 100 and np.sum(highSn) > 0:
                    ax.plot(
                        cast(Vector, xs[highSn]),
                        cast(Vector, ys[highSn]),
                        marker="x",
                        ms=4,
                        mec="w",
                        mew=2,
                        ls="none",
                    )
                    (highSnLine,) = ax.plot(
                        cast(Vector, xs[highSn]),
                        cast(Vector, ys[highSn]),
                        color=color,
                        marker="x",
                        ms=4,
                        ls="none",
                        label="High SN",
                    )
                    linesForLegend.append(highSnLine)
                    xMin = np.min(cast(Vector, xs[highSn]))
                else:
                    ax.axvline(highStats.approxMag, color=color, ls="--")

                if np.sum(lowSn) < 100 and np.sum(lowSn) > 0:
                    ax.plot(
                        cast(Vector, xs[lowSn]),
                        cast(Vector, ys[lowSn]),
                        marker="+",
                        ms=4,
                        mec="w",
                        mew=2,
                        ls="none",
                    )
                    (lowSnLine,) = ax.plot(
                        cast(Vector, xs[lowSn]),
                        cast(Vector, ys[lowSn]),
                        color=color,
                        marker="+",
                        ms=4,
                        ls="none",
                        label="Low SN",
                    )
                    linesForLegend.append(lowSnLine)
                    if xMin is None or xMin > np.min(cast(Vector, xs[lowSn])):
                        xMin = np.min(cast(Vector, xs[lowSn]))
                else:
                    ax.axvline(lowStats.approxMag, color=color, ls=":")

            else:
                ax.plot(xs, ys, ".", ms=5, alpha=0.3, mfc=color, mec=color, zorder=-1)
                meds = np.array([np.nanmedian(ys)] * len(xs))
                (medLine,) = ax.plot(xs, meds, color, label=f"Median: {np.nanmedian(ys):0.3g}", lw=0.8)
                linesForLegend.append(medLine)
                sigMads = np.array([sigmaMad(ys, nan_policy="omit")] * len(xs))
                (sigMadLine,) = ax.plot(
                    xs,
                    meds + 1.0 * sigMads,
                    color,
                    alpha=0.8,
                    lw=0.8,
                    label=r"$\sigma_{MAD}$: " + f"{sigMads[0]:0.3g}",
                )
                ax.plot(xs, meds - 1.0 * sigMads, color, alpha=0.8)
                linesForLegend.append(sigMadLine)
                histIm = None

        # Set the scatter plot limits
        # TODO: Make this not work by accident
        if len(cast(Vector, data["yStars"])) > 0:
            plotMed = np.nanmedian(cast(Vector, data["yStars"]))
        else:
            plotMed = np.nanmedian(cast(Vector, data["yGalaxies"]))
        # Ignore types below pending making this not working my accident
        if len(xs) < 2:  # type: ignore
            meds = [np.median(ys)]  # type: ignore
        if self.yLims:
            ax.set_ylim(self.yLims[0], self.yLims[1])  # type: ignore
        else:
            numSig = 4
            yLimMin = plotMed - numSig * sigMadYs  # type: ignore
            yLimMax = plotMed + numSig * sigMadYs  # type: ignore
            while (yLimMax < np.max(meds) or yLimMin > np.min(meds)) and numSig < 10:  # type: ignore
                numSig += 1

            numSig += 1
            yLimMin = plotMed - numSig * sigMadYs  # type: ignore
            yLimMax = plotMed + numSig * sigMadYs  # type: ignore
            ax.set_ylim(yLimMin, yLimMax)

        if self.xLims:
            ax.set_xlim(self.xLims[0], self.xLims[1])  # type: ignore
        elif len(xs) > 2:  # type: ignore
            if xMin is None:
                xMin = xs1 - 2 * xScale  # type: ignore
            ax.set_xlim(xMin, xs97 + 2 * xScale)  # type: ignore

        # Add a line legend
        ax.legend(
            handles=linesForLegend,
            ncol=4,
            fontsize=6,
            loc="upper left",
            framealpha=0.9,
            edgecolor="k",
            borderpad=0.4,
            handlelength=1,
        )

        # Add axes labels
        ax.set_ylabel(self.yAxisLabel, fontsize=10, labelpad=10)
        ax.set_xlabel(self.xAxisLabel, fontsize=10, labelpad=2)

        return ax, histIm

    def _makeTopHistogram(
        self, data: KeyedData, figure: Figure, gs: gridspec.GridSpec, ax: Axes, **kwargs
    ) -> None:
        # Top histogram
        totalX: list[Vector] = []
        if "stars" in self.plotTypes:  # type: ignore
            totalX.append(cast(Vector, data["xStars"]))
        if "galaxies" in self.plotTypes:  # type: ignore
            totalX.append(cast(Vector, data["xGalaxies"]))
        if "unknown" in self.plotTypes:  # type: ignore
            totalX.append(cast(Vector, data["xUknown"]))
        if "any" in self.plotTypes:  # type: ignore
            totalX.append(cast(Vector, data["x"]))

        totalXChained = [x for x in chain.from_iterable(totalX) if x == x]

        topHist = figure.add_subplot(gs[0, :-1], sharex=ax)
        topHist.hist(
            totalXChained, bins=100, color="grey", alpha=0.3, log=True, label=f"All ({len(totalXChained)})"
        )
        if "galaxies" in self.plotTypes:  # type: ignore
            topHist.hist(
                data["xGalaxies"],
                bins=100,
                color="firebrick",
                histtype="step",
                log=True,
                label=f"Galaxies ({len(cast(Vector, data['xGalaxies']))})",
            )
        if "stars" in self.plotTypes:  # type: ignore
            topHist.hist(
                data["xStars"],
                bins=100,
                color="midnightblue",
                histtype="step",
                log=True,
                label=f"Stars ({len(cast(Vector, data['xStars']))})",
            )
        topHist.axes.get_xaxis().set_visible(False)
        topHist.set_ylabel("Number", fontsize=8)
        topHist.legend(fontsize=6, framealpha=0.9, borderpad=0.4, loc="lower left", ncol=3, edgecolor="k")

        # Side histogram

    def _makeSideHistogram(
        self,
        data: KeyedData,
        figure: Figure,
        gs: gridspec.Gridspec,
        ax: Axes,
        histIm: Optional[PolyCollection],
        **kwargs,
    ) -> None:
        sideHist = figure.add_subplot(gs[1:, -1], sharey=ax)

        totalY: list[Vector] = []
        if "stars" in self.plotTypes:  # type: ignore
            totalY.append(cast(Vector, data["yStars"]))
        if "galaxies" in self.plotTypes:  # type: ignore
            totalY.append(cast(Vector, data["yGalaxies"]))
        if "unknown" in self.plotTypes:  # type: ignore
            totalY.append(cast(Vector, data["yUknown"]))
        if "any" in self.plotTypes:  # type: ignore
            totalY.append(cast(Vector, data["y"]))
        totalYChained = [y for y in chain.from_iterable(totalY) if y == y]

        # cheat to get the total count while iterating once
        yLimMin, yLimMax = ax.get_ylim()
        bins = np.linspace(yLimMin, yLimMax)
        sideHist.hist(
            totalYChained,
            bins=bins,
            color="grey",
            alpha=0.3,
            orientation="horizontal",
            log=True,
        )
        if "galaxies" in self.plotTypes:  # type: ignore
            sideHist.hist(
                [g for g in cast(Vector, data["yGalaxies"]) if g == g],
                bins=bins,
                color="firebrick",
                histtype="step",
                orientation="horizontal",
                log=True,
            )
            sideHist.hist(
                cast(Vector, data["yGalaxies"])[cast(Vector, data["galaxiesHighSNMask"])],
                bins=bins,
                color="firebrick",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls="--",
            )
            sideHist.hist(
                cast(Vector, data["yGalaxies"])[cast(Vector, data["galaxiesLowSNMask"])],
                bins=bins,
                color="firebrick",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls=":",
            )

        if "stars" in self.plotTypes:  # type: ignore
            sideHist.hist(
                [s for s in cast(Vector, data["yStars"]) if s == s],
                bins=bins,
                color="midnightblue",
                histtype="step",
                orientation="horizontal",
                log=True,
            )
            sideHist.hist(
                cast(Vector, data["yStars"])[cast(Vector, data["starsHighSNMask"])],
                bins=bins,
                color="midnightblue",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls="--",
            )
            sideHist.hist(
                cast(Vector, data["yStars"])[cast(Vector, data["starsLowSNMask"])],
                bins=bins,
                color="midnightblue",
                histtype="step",
                orientation="horizontal",
                log=True,
                ls=":",
            )

        sideHist.axes.get_yaxis().set_visible(False)
        sideHist.set_xlabel("Number", fontsize=8)
        if self.plot2DHist and histIm is not None:
            divider = make_axes_locatable(sideHist)
            cax = divider.append_axes("right", size="8%", pad=0)
            plt.colorbar(histIm, cax=cax, orientation="vertical", label="Number of Points Per Bin")
