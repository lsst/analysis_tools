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

__all__ = ("ScatterPlotStatsAction", "ScatterPlotWithTwoHists")

from itertools import chain
from typing import Mapping, NamedTuple, Optional, cast

import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField
from lsst.pex.config.listField import ListField
from lsst.skymap import BaseSkyMap
from matplotlib import gridspec
from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from matplotlib.figure import Figure
from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, PlotAction, Scalar, ScalarType, Vector
from ...statistics import nansigmaMad, sigmaMad
from ..keyedData.summaryStatistics import SummaryStatisticAction
from ..scalar import MedianAction
from ..vector import MagColumnNanoJansky, SnSelector
from .plotUtils import addPlotInfo, addSummaryPlot, generateSummaryStats, mkColormap

# ignore because coolwarm is actually part of module
cmapPatch = plt.cm.coolwarm.copy()  # type: ignore
cmapPatch.set_bad(color="none")


class ScatterPlotStatsAction(KeyedDataAction):
    """Calculates the statistics needed for the 
    scatter plot with two hists.
    """

    vectorKey = Field[str](doc="Vector on which to compute statistics")
    highSNSelector = ConfigurableActionField[SnSelector](
        doc="Selector used to determine high SN Objects", default=SnSelector(threshold=2700)
    )
    lowSNSelector = ConfigurableActionField[SnSelector](
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
        highMaskKey = f'{(self.identity or "").lower()}HighSNMask'
        results[highMaskKey] = self.highSNSelector(data, **kwargs)

        lowMaskKey = f'{(self.identity or "").lower()}LowSNMask'
        results[lowMaskKey] = self.lowSNSelector(data, **kwargs)

        prefix = f"{band}_" if (band := kwargs.get("band")) else ""
        fluxes = data[self.fluxType.format(band=band)] if band is not None else None

        statAction = SummaryStatisticAction(vectorKey=self.vectorKey)

        # this is sad, but pex_config seems to have broken behavior that
        # is dangerous to fix
        statAction.setDefaults()

        medianAction = MedianAction(vectorKey="mag")
        magAction = MagColumnNanoJansky(vectorKey="flux")

        for maskKey, binName in ((lowMaskKey, "low"), (highMaskKey, "high")):
            name = f"{prefix}{binName}SN{self.identity.capitalize() if self.identity else ''}"
            # set the approxMag to the median mag in the SN selection
            results[f"{name}_approxMag".format(**kwargs)] = (
                medianAction({"mag": magAction({"flux": fluxes[results[maskKey]]})})  # type: ignore
                if band is not None
                else np.nan
            )
            stats = statAction(data, **(kwargs | {"mask": results[maskKey]})).items()
            for suffix, value in stats:
                tmpKey = f"{name}_{suffix}".format(**kwargs)
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
    """Makes a scatter plot of the data with a marginal 
    histogram for each axis.
    """

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

    addSummaryPlot = Field[bool](
        doc="Add a summary plot to the figure?",
        default=False,
    )

    _stats = ("median", "sigmaMad", "count", "approxMag")

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | ScalarType]] = []
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
            base.append(("xUnknown", Vector))
            base.append(("yUnknown", Vector))
            base.append(("unknownHighSNMask", Vector))
            base.append(("unknownLowSNMask", Vector))
            # statistics
            for name in self._stats:
                base.append((f"{{band}}_highSNUnknown_{name}", Scalar))
                base.append((f"{{band}}_lowSNUnknown_{name}", Scalar))
        if "any" in self.plotTypes:  # type: ignore
            base.append(("x", Vector))
            base.append(("y", Vector))
            base.append(("anyHighSNMask", Vector))
            base.append(("anySNMask", Vector))
            # statistics
            for name in self._stats:
                base.append((f"{{band}}_highSNAny_{name}", Scalar))
                base.append((f"{{band}}_lowSNAny_{name}", Scalar))
        base.append(("lowSnThreshold", Scalar))
        base.append(("highSnThreshold", Scalar))

        if self.addSummaryPlot:
            base.append(("patch", Vector))

        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

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
        skymap: BaseSkyMap,
        plotInfo: Mapping[str, str],
        sumStats: Optional[Mapping] = None,
        **kwargs,
    ) -> Figure:
        """Makes a generic plot with a 2D histogram and collapsed histograms of
        each axis.

        Parameters
        ----------
        data : `KeyedData`
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

        If this function is being used within the pipetask framework 
        that takes care of making sure that data has all the required 
        elements but if you are runnign this as a standalone function 
        then you will need to provide the following things in the 
        input data.

        If stars is in self.plotTypes:
        xStars, yStars, starsHighSNMask, starsLowSNMask and 
        {band}_highSNStars_{name}, {band}_lowSNStars_{name} 
        where name is median, sigma_Mad, count and approxMag.

        If it is for galaxies/unknowns then replace stars in 
        the above names with galaxies/unknowns.

        if it is for any (which covers all the points) then it 
        becomes, x, y, and any instead of stars for the other 
        parameters given above.

        In every case it is expected that data contains:
        lowSnThreshold, highSnThreshold and patch 
        (if the summary plot is being plotted).
        """
        if not self.plotTypes:
            noDataFig = Figure()
            noDataFig.text(0.3, 0.5, "No data to plot after selectors applied")
            noDataFig = addPlotInfo(noDataFig, plotInfo)
            return noDataFig

        # Set default color and line style for the horizontal
        # reference line at 0
        if "hlineColor" not in kwargs:
            kwargs["hlineColor"] = "black"

        if "hlineStyle" not in kwargs:
            kwargs["hlineStyle"] = (0, (1, 4))

        fig = plt.figure(dpi=300)
        gs = gridspec.GridSpec(4, 4)

        # add the various plot elements
        ax, imhist = self._scatterPlot(data, fig, gs, **kwargs)
        self._makeTopHistogram(data, fig, gs, ax, **kwargs)
        self._makeSideHistogram(data, fig, gs, ax, imhist, **kwargs)
        # Needs info from run quantum
        if self.addSummaryPlot:
            sumStats = generateSummaryStats(data, skymap, plotInfo)
            label = self.yAxisLabel
            fig = addSummaryPlot(fig, gs[0, -1], sumStats, label)

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
        for j, (xs, ys, highSn, lowSn, color, cmap, highStats, lowStats) in enumerate(toPlotList):
            highSn = cast(Vector, highSn)
            lowSn = cast(Vector, lowSn)
            # ensure the columns are actually array
            xs = np.array(xs)
            ys = np.array(ys)
            sigMadYs = nansigmaMad(ys)
            if len(xs) < 2:
                (medLine,) = ax.plot(
                    xs, np.nanmedian(ys), color, label=f"Median: {np.nanmedian(ys):.2g}", lw=0.8
                )
                linesForLegend.append(medLine)
                sigMads = np.array([nansigmaMad(ys)] * len(xs))
                (sigMadLine,) = ax.plot(
                    xs,
                    np.nanmedian(ys) + 1.0 * sigMads,
                    color,
                    alpha=0.8,
                    lw=0.8,
                    label=r"$\sigma_{MAD}$: " + f"{sigMads[0]:.2g}",
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

                for i, xEdge in enumerate(xEdgesPlot):
                    ids = np.where((xs < xEdge) & (xs > xEdges[i]) & (np.isfinite(ys)))[0]
                    med = np.nanmedian(ys[ids])
                    sigMad = sigmaMad(ys[ids], nan_policy="omit")
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
                statText = f"S/N > {highThresh:0.4g} Stats ({self.magLabel} < {highStats.approxMag:0.4g})\n"
                highStatsStr = (
                    f"Median: {highStats.median:0.4g}    "
                    + r"$\sigma_{MAD}$: "
                    + f"{highStats.sigmaMad:0.4g}    "
                    + r"N$_{points}$: "
                    + f"{highStats.count}"
                )
                statText += highStatsStr
                fig.text(xPos, 0.090, statText, bbox=bbox, transform=fig.transFigure, fontsize=6)

                bbox = dict(edgecolor=color, linestyle=":", facecolor="none")
                lowThresh = data["lowSnThreshold"]
                statText = f"S/N > {lowThresh:0.4g} Stats ({self.magLabel} < {lowStats.approxMag:0.4g})\n"
                lowStatsStr = (
                    f"Median: {lowStats.median:0.4g}    "
                    + r"$\sigma_{MAD}$: "
                    + f"{lowStats.sigmaMad:0.4g}    "
                    + r"N$_{points}$: "
                    + f"{lowStats.count}"
                )
                statText += lowStatsStr
                fig.text(xPos, 0.020, statText, bbox=bbox, transform=fig.transFigure, fontsize=6)

                if self.plot2DHist:
                    histIm = ax.hexbin(xs[inside], ys[inside], gridsize=75, cmap=cmap, mincnt=1, zorder=-3)

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
                sigMads = np.array([nansigmaMad(ys)] * len(xs))
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

        # Add a horizontal reference line at 0 to the scatter plot
        ax.axhline(0, color=kwargs["hlineColor"], ls=kwargs["hlineStyle"], alpha=0.7, zorder=-2)

        # Set the scatter plot limits
        # TODO: Make this not work by accident
        if len(cast(Vector, data["yStars"])) > 0:
            plotMed = np.nanmedian(cast(Vector, data["yStars"]))
        else:
            plotMed = np.nanmedian(cast(Vector, data["yGalaxies"]))
        # Ignore types below pending making this not working my accident
        if len(xs) < 2:  # type: ignore
            meds = [np.nanmedian(ys)]  # type: ignore
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

        # Add a horizontal reference line at 0 to the side histogram
        sideHist.axhline(0, color=kwargs["hlineColor"], ls=kwargs["hlineStyle"], alpha=0.7, zorder=-2)

        sideHist.axes.get_yaxis().set_visible(False)
        sideHist.set_xlabel("Number", fontsize=8)
        if self.plot2DHist and histIm is not None:
            divider = make_axes_locatable(sideHist)
            cax = divider.append_axes("right", size="8%", pad=0)
            plt.colorbar(histIm, cax=cax, orientation="vertical", label="Number of Points Per Bin")
