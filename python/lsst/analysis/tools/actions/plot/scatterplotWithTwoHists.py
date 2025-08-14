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

import math
from typing import Mapping, NamedTuple, Optional, cast

import matplotlib.colors
import matplotlib.patheffects as pathEffects
import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField
from lsst.pex.config.listField import ListField
from lsst.utils.plotting import (
    galaxies_cmap,
    galaxies_color,
    make_figure,
    set_rubin_plotstyle,
    stars_cmap,
    stars_color,
)
from matplotlib import gridspec
from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from matplotlib.figure import Figure
from matplotlib.path import Path
from matplotlib.ticker import LogFormatterMathtext, NullFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, PlotAction, Scalar, ScalarType, Vector
from ...math import nanMedian, nanSigmaMad
from ..keyedData.summaryStatistics import SummaryStatisticAction
from ..scalar import MedianAction
from ..vector import ConvertFluxToMag, SnSelector
from .plotUtils import addPlotInfo, addSummaryPlot, generateSummaryStats


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
    prefix = Field[str](
        doc="Prefix for all output fields; will use self.identity if None",
        optional=True,
        default=None,
    )
    fluxType = Field[str](doc="Vector key to use to compute signal to noise ratio", default="{band}_psfFlux")
    suffix = Field[str](doc="Suffix for all output fields", default="")

    def _get_key_prefix(self):
        prefix = self.prefix if self.prefix else (self.identity if self.identity else "")
        return prefix

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield (self.vectorKey, Vector)
        yield (self.fluxType, Vector)
        yield from self.highSNSelector.getInputSchema()
        yield from self.lowSNSelector.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        prefix = self._get_key_prefix()
        prefix_lower = prefix.lower() if prefix else ""
        prefix_upper = prefix.capitalize() if prefix else ""
        suffix = self.suffix
        return (
            (f"{prefix_lower}HighSNMask{suffix}", Vector),
            (f"{prefix_lower}LowSNMask{suffix}", Vector),
            (f"{{band}}_lowSN{prefix_upper}_median{suffix}", Scalar),
            (f"{{band}}_lowSN{prefix_upper}_sigmaMad{suffix}", Scalar),
            (f"{{band}}_lowSN{prefix_upper}_count{suffix}", Scalar),
            (f"{{band}}_lowSN{prefix_upper}_approxMag{suffix}", Scalar),
            (f"{{band}}_highSN{prefix_upper}_median{suffix}", Scalar),
            (f"{{band}}_highSN{prefix_upper}_sigmaMad{suffix}", Scalar),
            (f"{{band}}_highSN{prefix_upper}_count{suffix}", Scalar),
            (f"{{band}}_highSN{prefix_upper}_approxMag{suffix}", Scalar),
            (f"{prefix_lower}HighSNThreshold{suffix}", Scalar),
            (f"{prefix_lower}LowSNThreshold{suffix}", Scalar),
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        prefix = self._get_key_prefix()
        prefix_lower = prefix.lower() if prefix else ""
        prefix_upper = prefix.capitalize() if prefix else ""
        suffix = self.suffix
        highMaskKey = f"{prefix_lower}HighSNMask{suffix}"
        results[highMaskKey] = self.highSNSelector(data, **kwargs)

        lowMaskKey = f"{prefix_lower}LowSNMask{suffix}"
        results[lowMaskKey] = self.lowSNSelector(data, **kwargs)

        prefix_band = f"{band}_" if (band := kwargs.get("band")) else ""
        fluxes = data[self.fluxType.format(band=band)] if band is not None else None

        statAction = SummaryStatisticAction(vectorKey=self.vectorKey)

        # this is sad, but pex_config seems to have broken behavior that
        # is dangerous to fix
        statAction.setDefaults()

        medianAction = MedianAction(vectorKey="mag")
        magAction = ConvertFluxToMag(vectorKey="flux")

        for maskKey, binName in ((lowMaskKey, "low"), (highMaskKey, "high")):
            name = f"{prefix_band}{binName}SN{prefix_upper}"
            # set the approxMag to the median mag in the SN selection
            results[f"{name}_approxMag{suffix}".format(**kwargs)] = (
                medianAction({"mag": magAction({"flux": fluxes[results[maskKey]]})})  # type: ignore
                if band is not None
                else np.nan
            )
            stats = statAction(data, **(kwargs | {"mask": results[maskKey]})).items()
            for name_stat, value in stats:
                tmpKey = f"{name}_{name_stat}{suffix}".format(**kwargs)
                results[tmpKey] = value
        results[f"{prefix_lower}HighSNThreshold{suffix}"] = self.highSNSelector.threshold  # type: ignore
        results[f"{prefix_lower}LowSNThreshold{suffix}"] = self.lowSNSelector.threshold  # type: ignore

        return results


def _validObjectTypes(value):
    return value in ("stars", "galaxies", "unknown", "any")


# ignore type because of conflicting name on tuple baseclass
class _StatsContainer(NamedTuple):
    median: Scalar
    sigmaMad: Scalar
    count: Scalar  # type: ignore
    approxMag: Scalar


class DataTypeDefaults(NamedTuple):
    suffix_stat: str
    suffix_xy: str
    color: str
    colormap: matplotlib.colors.Colormap | None


class LogFormatterExponentSci(LogFormatterMathtext):
    """
    Format values following scientific notation.

    Unlike the matplotlib LogFormatterExponent, this will print near-integer
    coefficients with a base between 0 and 2 such as 500 as 500 (if base10)
    or 5e2 otherwise.
    """

    def _non_decade_format(self, sign_string, base, fx, usetex):
        """Return string for non-decade locations."""
        b = float(base)
        exponent = math.floor(fx)
        coeff = b ** (fx - exponent)
        rounded = round(coeff)
        if math.isclose(coeff, rounded):
            if (base == "10") and (0 <= exponent <= 3):
                return f"{sign_string}{rounded}{'0'*int(exponent)}"
            coeff = rounded
        return f"{sign_string}{coeff:1.1f}e{exponent}"


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

    legendLocation = Field[str](doc="Legend position within main plot", default="upper left")
    nBins = Field[float](doc="Number of bins on x axis", default=40.0)
    plot2DHist = Field[bool](
        doc="Plot a 2D histogram in dense areas of points on the scatter plot."
        "Doesn't look great if plotting multiple datasets on top of each other.",
        default=True,
    )
    plotTypes = ListField[str](
        doc="Selection of types of objects to plot. Can take any combination of"
        " stars, galaxies, unknown, any",
        optional=False,
        itemCheck=_validObjectTypes,
    )

    addSummaryPlot = Field[bool](
        doc="Add a summary plot to the figure?",
        default=True,
    )
    histMinimum = Field[float](
        doc="Minimum value for the histogram count axis",
        default=0.3,
    )
    xHistMaxLabels = Field[int](
        doc="Maximum number of labels for ticks on the x-axis marginal histogram",
        default=3,
        check=lambda x: x >= 2,
    )
    yHistMaxLabels = Field[int](
        doc="Maximum number of labels for ticks on the y-axis marginal histogram",
        default=3,
        check=lambda x: x >= 2,
    )

    suffix_x = Field[str](doc="Suffix for all x-axis action inputs", optional=True, default="")
    suffix_y = Field[str](doc="Suffix for all y-axis action inputs", optional=True, default="")
    suffix_stat = Field[str](doc="Suffix for all binned statistic action inputs", optional=True, default="")

    publicationStyle = Field[bool](doc="Slimmed down publication style plot?", default=False)

    _stats = ("median", "sigmaMad", "count", "approxMag")
    _datatypes = {
        "galaxies": DataTypeDefaults(
            suffix_stat="Galaxies",
            suffix_xy="Galaxies",
            color=galaxies_color(),
            colormap=galaxies_cmap(single_color=True),
        ),
        "stars": DataTypeDefaults(
            suffix_stat="Stars",
            suffix_xy="Stars",
            color=stars_color(),
            colormap=stars_cmap(single_color=True),
        ),
        "unknown": DataTypeDefaults(
            suffix_stat="Unknown",
            suffix_xy="Unknown",
            color="green",
            colormap=None,
        ),
        "any": DataTypeDefaults(
            suffix_stat="Any",
            suffix_xy="",
            color="purple",
            colormap=None,
        ),
    }

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | ScalarType]] = []
        for name_datatype in self.plotTypes:
            config_datatype = self._datatypes[name_datatype]
            if not self.publicationStyle:
                base.append((f"x{config_datatype.suffix_xy}{self.suffix_x}", Vector))
                base.append((f"y{config_datatype.suffix_xy}{self.suffix_y}", Vector))
                base.append((f"{name_datatype}HighSNMask{self.suffix_stat}", Vector))
                base.append((f"{name_datatype}LowSNMask{self.suffix_stat}", Vector))
                # statistics
                for name in self._stats:
                    base.append(
                        (f"{{band}}_highSN{config_datatype.suffix_stat}_{name}{self.suffix_stat}", Scalar)
                    )
                    base.append(
                        (f"{{band}}_lowSN{config_datatype.suffix_stat}_{name}{self.suffix_stat}", Scalar)
                    )
                base.append((f"{name_datatype}LowSNThreshold{self.suffix_stat}", Scalar))
                base.append((f"{name_datatype}HighSNThreshold{self.suffix_stat}", Scalar))

        if self.addSummaryPlot and not self.publicationStyle:
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
            raise ValueError(
                f"Task needs keys {remainder} but they were not found in input keys" f" {list(data.keys())}"
            )
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
        """Makes a generic plot with a 2D histogram and collapsed histograms of
        each axis.

        Parameters
        ----------
        data : `KeyedData`
            The catalog to plot the points from.
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:

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
        elements but if you are running this as a standalone function
        then you will need to provide the following things in the
        input data.

        * If stars is in self.plotTypes:
            xStars, yStars, starsHighSNMask, starsLowSNMask and
            {band}_highSNStars_{name}, {band}_lowSNStars_{name}
            where name is median, sigma_Mad, count and approxMag.

        * If it is for galaxies/unknowns then replace stars in the above
          names with galaxies/unknowns.

        * If it is for any (which covers all the points) then it
          becomes, x, y, and any instead of stars for the other
          parameters given above.

        * In every case it is expected that data contains:
            lowSnThreshold, highSnThreshold and patch
            (if the summary plot is being plotted).

        Examples
        --------
        An example of the plot produced from this code is here:

        .. image:: /_static/analysis_tools/scatterPlotExample.png

        For a detailed example of how to make a plot from the command line
        please see the
        :ref:`getting started guide<analysis-tools-getting-started>`.
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

        set_rubin_plotstyle()
        fig = make_figure()
        gs = gridspec.GridSpec(4, 4)

        # add the various plot elements
        ax, imhist = self._scatterPlot(data, fig, gs, **kwargs)
        if ax is None:
            noDataFig = Figure()
            noDataFig.text(0.3, 0.5, "No data to plot after selectors applied")
            if not self.publicationStyle:
                noDataFig = addPlotInfo(noDataFig, plotInfo)
            return noDataFig

        self._makeTopHistogram(data, fig, gs, ax, **kwargs)
        self._makeSideHistogram(data, fig, gs, ax, imhist, **kwargs)
        # Needs info from run quantum
        skymap = kwargs.get("skymap", None)
        if self.addSummaryPlot and skymap is not None and not self.publicationStyle:
            sumStats = generateSummaryStats(data, skymap, plotInfo)
            label = self.yAxisLabel
            fig = addSummaryPlot(fig, gs[0, -1], sumStats, label)

        fig.canvas.draw()
        # TODO: Check if these spacings can be defined less arbitrarily
        fig.subplots_adjust(
            wspace=0.0,
            hspace=0.0,
            bottom=0.13 if self.publicationStyle else 0.22,
            left=0.18 if self.publicationStyle else 0.21,
            right=0.92 if self.publicationStyle else None,
            top=0.98 if self.publicationStyle else None,
        )
        if not self.publicationStyle:
            fig = addPlotInfo(fig, plotInfo)
        return fig

    def _scatterPlot(
        self, data: KeyedData, fig: Figure, gs: gridspec.GridSpec, **kwargs
    ) -> tuple[Axes, Optional[PolyCollection]]:
        suf_x = self.suffix_x
        suf_y = self.suffix_y
        suf_stat = self.suffix_stat
        # Main scatter plot
        ax = fig.add_subplot(gs[1:, :-1])

        binThresh = 5

        yBinsOut = []
        linesForLegend = []

        toPlotList = []
        histIm = None
        highStats: _StatsContainer
        lowStats: _StatsContainer

        for name_datatype in self.plotTypes:
            config_datatype = self._datatypes[name_datatype]
            highArgs = {}
            lowArgs = {}
            if not self.publicationStyle:
                for name in self._stats:
                    highArgs[name] = cast(
                        Scalar,
                        data[
                            f"{{band}}_highSN{config_datatype.suffix_stat}_{name}{suf_stat}".format(**kwargs)
                        ],
                    )
                    lowArgs[name] = cast(
                        Scalar,
                        data[
                            f"{{band}}_lowSN{config_datatype.suffix_stat}_{name}{suf_stat}".format(**kwargs)
                        ],
                    )
                highStats = _StatsContainer(**highArgs)
                lowStats = _StatsContainer(**lowArgs)

                toPlotList.append(
                    (
                        data[f"x{config_datatype.suffix_xy}{suf_x}"],
                        data[f"y{config_datatype.suffix_xy}{suf_y}"],
                        data[f"{name_datatype}HighSNMask{suf_stat}"],
                        data[f"{name_datatype}LowSNMask{suf_stat}"],
                        data[f"{name_datatype}HighSNThreshold{suf_stat}"],
                        data[f"{name_datatype}LowSNThreshold{suf_stat}"],
                        config_datatype.color,
                        config_datatype.colormap,
                        highStats,
                        lowStats,
                    )
                )
            else:
                toPlotList.append(
                    (
                        data[f"x{config_datatype.suffix_xy}{suf_x}"],
                        data[f"y{config_datatype.suffix_xy}{suf_y}"],
                        [],
                        [],
                        [],
                        [],
                        config_datatype.color,
                        config_datatype.colormap,
                        [],
                        [],
                    )
                )

        xLims = self.xLims if self.xLims is not None else [np.inf, -np.inf]

        # If there is no data to plot make a
        # no data figure
        numData = 0
        for xs, _, _, _, _, _, _, _, _, _ in toPlotList:
            numData += len(xs)
        if numData == 0:
            return None, None

        for j, (
            xs,
            ys,
            highSn,
            lowSn,
            highThresh,
            lowThresh,
            color,
            cmap,
            highStats,
            lowStats,
        ) in enumerate(toPlotList):
            highSn = cast(Vector, highSn)
            lowSn = cast(Vector, lowSn)
            # ensure the columns are actually array
            xs = np.array(xs)
            ys = np.array(ys)
            sigMadYs = nanSigmaMad(ys)
            # plot lone median point if there's not enough data to measure more
            n_xs = len(xs)
            if n_xs == 0 or not np.isfinite(sigMadYs):
                continue
            elif n_xs < 10:
                xs = [nanMedian(xs)]
                sigMads = np.array([nanSigmaMad(ys)])
                ys = np.array([nanMedian(ys)])
                (medLine,) = ax.plot(xs, ys, color, label=f"Median: {ys[0]:.2g}", lw=0.8)
                linesForLegend.append(medLine)
                (sigMadLine,) = ax.plot(
                    xs,
                    ys + 1.0 * sigMads,
                    color,
                    alpha=0.8,
                    lw=0.8,
                    label=r"$\sigma_{MAD}$: " + f"{sigMads[0]:.2g}",
                )
                ax.plot(xs, ys - 1.0 * sigMads, color, alpha=0.8)
                linesForLegend.append(sigMadLine)
                histIm = None
                continue

            if self.xLims:
                xMin, xMax = self.xLims
            else:
                # Chop off 1/3% from only the finite xs values
                # (there may be +/-np.Inf values)
                # TODO: This should be configurable
                # but is there a good way to avoid redundant config params
                # without using slightly annoying subconfigs?
                xs1, xs97 = np.nanpercentile(xs[np.isfinite(xs)], (1, 97))
                xScale = (xs97 - xs1) / 20.0  # This is ~5% of the data range
                xMin, xMax = (xs1 - xScale, xs97 + xScale)
                xLims[0] = min(xLims[0], xMin)
                xLims[1] = max(xLims[1], xMax)

            xEdges = np.arange(xMin, xMax, (xMax - xMin) / self.nBins)
            medYs = nanMedian(ys)
            fiveSigmaHigh = medYs + 5.0 * sigMadYs
            fiveSigmaLow = medYs - 5.0 * sigMadYs
            binSize = (fiveSigmaHigh - fiveSigmaLow) / 101.0
            # When the binsize is 0 try using the 1st and 99th
            # percentile instead of the sigmas.
            if binSize == 0.0:
                p1, p99 = np.nanpercentile(ys, [1, 99])
                binSize = (p99 - p1) / 101.0

            # If fiveSigmaHigh and fiveSigmaLow are the same
            # then use the 1st and 99th percentiles to define the
            # yEdges.
            yEdges = np.arange(fiveSigmaLow, fiveSigmaHigh, binSize)
            if fiveSigmaLow == fiveSigmaHigh:
                yEdges = np.arange(p1, p99, binSize)

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
                    med = nanMedian(ys[ids])
                    sigMad = nanSigmaMad(ys[ids])
                    meds[i] = med
                    sigMads[i] = sigMad
                    threeSigMadVerts[i, :] = [xEdge, med + 3 * sigMad]
                    threeSigMadVerts[-(i + 1), :] = [xEdge, med - 3 * sigMad]

                if self.publicationStyle:
                    linecolor = "k"
                else:
                    linecolor = color

                (medLine,) = ax.plot(xEdgesPlot, meds, linecolor, label="Running Median")
                linesForLegend.append(medLine)

                # Make path to check which points lie within one sigma mad
                threeSigMadPath = Path(threeSigMadVerts, codes)

                if not self.publicationStyle:
                    # Add lines for the median +/- 3 * sigma MAD
                    (threeSigMadLine,) = ax.plot(
                        xEdgesPlot,
                        threeSigMadVerts[: len(xEdgesPlot), 1],
                        linecolor,
                        alpha=0.4,
                        label=r"3$\sigma_{MAD}$",
                    )
                    ax.plot(xEdgesPlot[::-1], threeSigMadVerts[len(xEdgesPlot) :, 1], linecolor, alpha=0.4)

                # Add lines for the median +/- 1 * sigma MAD
                (sigMadLine,) = ax.plot(
                    xEdgesPlot,
                    meds + 1.0 * sigMads,
                    linecolor,
                    alpha=0.8,
                    label=r"$\sigma_{MAD}$",
                    ls="dashed",
                )
                linesForLegend.append(sigMadLine)
                ax.plot(xEdgesPlot, meds - 1.0 * sigMads, linecolor, alpha=0.8, ls="dashed")

                if not self.publicationStyle:
                    # Add lines for the median +/- 2 * sigma MAD
                    (twoSigMadLine,) = ax.plot(
                        xEdgesPlot, meds + 2.0 * sigMads, linecolor, alpha=0.6, label=r"2$\sigma_{MAD}$"
                    )
                    linesForLegend.append(twoSigMadLine)
                    linesForLegend.append(threeSigMadLine)
                    ax.plot(xEdgesPlot, meds - 2.0 * sigMads, linecolor, alpha=0.6)

                # Check which points are outside 3 sigma MAD of the median
                # and plot these as points.
                inside = threeSigMadPath.contains_points(np.array([xs, ys]).T)
                ax.plot(xs[~inside], ys[~inside], ".", ms=5, alpha=0.2, mfc=color, mec="none", zorder=-1)

                if not self.publicationStyle:
                    # Add some stats text
                    xPos = 0.65 - 0.4 * j
                    bbox = dict(edgecolor=color, linestyle="--", facecolor="none")
                    statText = (
                        f"S/N > {highThresh:0.4g} Stats ({self.magLabel} < {highStats.approxMag:0.4g})\n"
                    )
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
                    extent = [xLims[0], xLims[1], self.yLims[0], self.yLims[1]] if self.yLims else None
                    histIm = ax.hexbin(
                        xs[inside],
                        ys[inside],
                        gridsize=75,
                        extent=extent,
                        cmap=cmap,
                        mincnt=1,
                        zorder=-3,
                        edgecolors=None,
                    )
                else:
                    ax.plot(xs[inside], ys[inside], ".", ms=3, alpha=0.2, mfc=color, mec=color, zorder=-1)

                if not self.publicationStyle:
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
                    else:
                        ax.axvline(lowStats.approxMag, color=color, ls=":")

            else:
                ax.plot(xs, ys, ".", ms=5, alpha=0.3, mfc=color, mec=color, zorder=-1)
                meds = np.array([nanMedian(ys)] * len(xs))
                (medLine,) = ax.plot(xs, meds, color, label=f"Median: {nanMedian(ys):0.3g}", lw=0.8)
                linesForLegend.append(medLine)
                sigMads = np.array([nanSigmaMad(ys)] * len(xs))
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
        suf_x = self.suffix_y
        # TODO: Make this not work by accident
        if f"yStars{suf_x}" in data and (len(cast(Vector, data[f"yStars{suf_x}"])) > 0):
            plotMed = nanMedian(cast(Vector, data[f"yStars{suf_x}"]))
        elif f"yGalaxies{suf_x}" in data and (len(cast(Vector, data[f"yGalaxies{suf_x}"])) > 0):
            plotMed = nanMedian(cast(Vector, data[f"yGalaxies{suf_x}"]))
        else:
            plotMed = np.nan

        # Ignore types below pending making this not working my accident
        if len(xs) < 2:  # type: ignore
            meds = [nanMedian(ys)]  # type: ignore
        if self.yLims:
            ax.set_ylim(self.yLims[0], self.yLims[1])  # type: ignore
        elif np.isfinite(plotMed):
            numSig = 4
            yLimMin = plotMed - numSig * sigMadYs  # type: ignore
            yLimMax = plotMed + numSig * sigMadYs  # type: ignore
            while (yLimMax < np.max(meds) or yLimMin > np.min(meds)) and numSig < 10:  # type: ignore
                numSig += 1

            numSig += 1
            yLimMin = plotMed - numSig * sigMadYs  # type: ignore
            yLimMax = plotMed + numSig * sigMadYs  # type: ignore
            ax.set_ylim(yLimMin, yLimMax)

        # This could be false if len(x) == 0 for xs in toPlotList
        # ... in which case nothing was plotted and limits are irrelevant
        if all(np.isfinite(xLims)):
            ax.set_xlim(xLims)

        # Add a line legend
        ax.legend(
            handles=linesForLegend,
            ncol=4,
            fontsize=6,
            loc=self.legendLocation,
            framealpha=0.9,
            edgecolor="k",
            borderpad=0.4,
            handlelength=3,
        )

        # Add axes labels
        band = kwargs.get("band", "unspecified")
        xlabel = self.xAxisLabel
        ylabel = self.yAxisLabel
        if "{band}" in xlabel:
            xlabel = xlabel.format(band=band)
        if "{band}" in ylabel:
            ylabel = ylabel.format(band=band)
        if self.publicationStyle:
            ax.set_ylabel(ylabel, labelpad=10)
            ax.set_xlabel(xlabel, labelpad=2)
        else:
            ax.set_ylabel(ylabel, labelpad=10, fontsize=10)
            ax.set_xlabel(xlabel, labelpad=2, fontsize=10)
            ax.tick_params(labelsize=8)

        return ax, histIm

    def _makeTopHistogram(
        self, data: KeyedData, figure: Figure, gs: gridspec.GridSpec, ax: Axes, **kwargs
    ) -> None:
        suf_x = self.suffix_x
        # Top histogram
        topHist = figure.add_subplot(gs[0, :-1], sharex=ax)
        x_min, x_max = ax.get_xlim()
        bins = np.linspace(x_min, x_max, 100)

        if "any" in self.plotTypes:
            x_any = f"x{self._datatypes['any'].suffix_xy}{suf_x}"
            keys_notany = [x for x in self.plotTypes if x != "any"]
        else:
            x_any = (
                np.concatenate([data[f"x{self._datatypes[key].suffix_xy}{suf_x}"] for key in self.plotTypes])
                if (len(self.plotTypes) > 1)
                else None
            )
            keys_notany = self.plotTypes
        if x_any is not None:
            if np.sum(x_any > 0) > 0:
                log = True
            else:
                log = False
            topHist.hist(x_any, bins=bins, color="grey", alpha=0.3, log=log, label=f"Any ({len(x_any)})")

        for key in keys_notany:
            config_datatype = self._datatypes[key]
            vector = np.array(data[f"x{config_datatype.suffix_xy}{suf_x}"])
            if np.sum(vector > 0) > 0:
                log = True
            else:
                log = False
            topHist.hist(
                vector,
                bins=bins,
                color=config_datatype.color,
                histtype="step",
                log=log,
                label=f"{config_datatype.suffix_stat} ({len(vector)})",
            )
        topHist.axes.get_xaxis().set_visible(False)
        topHist.set_ylabel("Count", fontsize=10 + 4 * self.publicationStyle)
        if not self.publicationStyle:
            topHist.legend(fontsize=6, framealpha=0.9, borderpad=0.4, loc="lower left", ncol=3, edgecolor="k")
            topHist.tick_params(labelsize=8)

        self._modifyHistogramTicks(topHist, do_x=False, max_labels=self.xHistMaxLabels)

    def _makeSideHistogram(
        self,
        data: KeyedData,
        figure: Figure,
        gs: gridspec.Gridspec,
        ax: Axes,
        histIm: Optional[PolyCollection],
        **kwargs,
    ) -> None:
        suf_y = self.suffix_y
        # Side histogram
        sideHist = figure.add_subplot(gs[1:, -1], sharey=ax)
        y_min, y_max = ax.get_ylim()
        bins = np.linspace(y_min, y_max, 100)

        if "any" in self.plotTypes:
            y_any = np.array(data[f"y{self._datatypes['any'].suffix_xy}{suf_y}"])
            keys_notany = [x for x in self.plotTypes if x != "any"]
        else:
            y_any = (
                np.concatenate(
                    [np.array(data[f"y{self._datatypes[key].suffix_xy}{suf_y}"]) for key in self.plotTypes]
                )
                if (len(self.plotTypes) > 1)
                else None
            )
            keys_notany = self.plotTypes
        if y_any is not None:
            sideHist.hist(
                np.array(y_any),
                bins=bins,
                color="grey",
                alpha=0.3,
                orientation="horizontal",
                log=np.any(y_any > 0),
            )

        kwargs_hist = dict(
            bins=bins,
            histtype="step",
            log=True,
            orientation="horizontal",
        )
        for key in keys_notany:
            config_datatype = self._datatypes[key]
            # If the data has no positive values then it
            # cannot be log scaled and it prints a bunch
            # of irritating warnings, in this case don't
            # try.
            numPos = np.sum(data[f"y{config_datatype.suffix_xy}{suf_y}"] > 0)

            if numPos <= 0:
                kwargs_hist["log"] = False

            vector = data[f"y{config_datatype.suffix_xy}{suf_y}"]
            sideHist.hist(
                vector,
                color=config_datatype.color,
                **kwargs_hist,
            )
            if not self.publicationStyle:
                sideHist.hist(
                    vector[cast(Vector, data[f"{key}HighSNMask{self.suffix_stat}"])],
                    color=config_datatype.color,
                    linestyle="--",
                    **kwargs_hist,
                )
                sideHist.hist(
                    vector[cast(Vector, data[f"{key}LowSNMask{self.suffix_stat}"])],
                    color=config_datatype.color,
                    **kwargs_hist,
                    linestyle=":",
                )

        # Add a horizontal reference line at 0 to the side histogram
        sideHist.axhline(0, color=kwargs["hlineColor"], ls=kwargs["hlineStyle"], alpha=0.7, zorder=-2)

        sideHist.axes.get_yaxis().set_visible(False)
        sideHist.set_xlabel("Count", fontsize=10 + 4 * self.publicationStyle)
        self._modifyHistogramTicks(sideHist, do_x=True, max_labels=self.yHistMaxLabels)

        if not self.publicationStyle:
            sideHist.tick_params(labelsize=8)
        if self.plot2DHist and histIm is not None:
            divider = make_axes_locatable(sideHist)
            cax = divider.append_axes("right", size="25%", pad=0)
            sideHist.get_figure().colorbar(histIm, cax=cax, orientation="vertical")
            text = cax.text(
                0.5,
                0.5,
                "Points Per Bin",
                color="k",
                rotation="vertical",
                transform=cax.transAxes,
                ha="center",
                va="center",
                fontsize=10,
            )
            text.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])

    def _modifyHistogramTicks(self, histogram, do_x: bool, max_labels: int):
        axis = histogram.get_xaxis() if do_x else histogram.get_yaxis()
        limits = list(histogram.get_xlim() if do_x else histogram.get_ylim())
        get_ticks = histogram.get_xticks if do_x else histogram.get_yticks
        ticks = get_ticks()
        # Let the minimum be larger then specified if the histogram has large
        # values everywhere, but cut it down a little so the lowest-valued bin
        # is still easily visible
        limits[0] = max(self.histMinimum, 0.9 * limits[0])
        # Round the upper limit to the nearest power of 10
        limits[1] = 10 ** (np.ceil(np.log10(limits[1]))) if (limits[1] > 0) else limits[1]
        for minor in (False, True):
            # Ignore ticks that are below the minimum value
            valid = (ticks >= limits[0]) & (ticks <= limits[1])
            labels = [label for label, _valid in zip(axis.get_ticklabels(minor=minor), valid) if _valid]
            if (n_labels := len(labels)) > max_labels:
                labels_new = [""] * n_labels
                # Skip the first label if we're not using minor axis labels
                # This helps avoid overlap with the scatter plot labels
                for idx_fill in np.round(np.linspace(1 - minor, n_labels - 1, max_labels)).astype(int):
                    labels_new[idx_fill] = labels[idx_fill]
                axis.set_ticks(ticks[valid], labels_new)
            # If there are enough major tick labels, disable minor tick labels
            if len(labels) >= 2:
                axis.set_minor_formatter(NullFormatter())
                break
            else:
                axis.set_minor_formatter(
                    LogFormatterExponentSci(minor_thresholds=(1, self.histMinimum / 10.0))
                )
                ticks = get_ticks(minor=True)

        (histogram.set_xlim if do_x else histogram.set_ylim)(limits)
