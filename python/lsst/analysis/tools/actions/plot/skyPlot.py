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

__all__ = ("SkyPlot",)

from typing import Mapping, Optional

import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import Field, ListField
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from ...statistics import nansigmaMad
from .plotUtils import addPlotInfo, mkColormap, plotProjectionWithBinning, sortAllArrays


class SkyPlot(PlotAction):
    """Plots the on sky distribution of a parameter.

    Plots the values of the parameter given for the z axis
    according to the positions given for x and y. Optimised
    for use with RA and Dec. Also calculates some basic
    statistics and includes those on the plot.
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis.", optional=False)
    yAxisLabel = Field[str](doc="Label to use for the y axis.", optional=False)
    zAxisLabel = Field[str](doc="Label to use for the z axis.", optional=False)

    plotOutlines = Field[bool](
        doc="Plot the outlines of the ccds/patches?",
        default=True,
    )

    plotTypes = ListField[str](
        doc="Selection of types of objects to plot. Can take any combination of"
        " stars, galaxies, unknown, mag, any.",
        optional=False,
        # itemCheck=_validatePlotTypes,
    )

    plotName = Field[str](doc="The name for the plot.", optional=False)

    fixAroundZero = Field[bool](
        doc="Fix the colorbar to be symmetric around zero.",
        default=False,
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        if "stars" in self.plotTypes:  # type: ignore
            base.append(("xStars", Vector))
            base.append(("yStars", Vector))
            base.append(("zStars", Vector))
            base.append(("starStatMask", Vector))
        if "galaxies" in self.plotTypes:  # type: ignore
            base.append(("xGalaxies", Vector))
            base.append(("yGalaxies", Vector))
            base.append(("zGalaxies", Vector))
            base.append(("galaxyStatMask", Vector))
        if "unknown" in self.plotTypes:  # type: ignore
            base.append(("xUnknowns", Vector))
            base.append(("yUnknowns", Vector))
            base.append(("zUnknowns", Vector))
            base.append(("unknownStatMask", Vector))
        if "any" in self.plotTypes:  # type: ignore
            base.append(("x", Vector))
            base.append(("y", Vector))
            base.append(("z", Vector))
            base.append(("statMask", Vector))

        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)
        # table is a dict that needs: x, y, run, skymap, filter, tract,

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """NOTE currently can only check that something is not a Scalar, not
        check that the data is consistent with Vector
        """
        needed = self.getInputSchema(**kwargs)
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in input")
        for name, typ in needed:
            isScalar = issubclass((colType := type(data[name.format(**kwargs)])), Scalar)
            if isScalar and typ != Scalar:
                raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}")

    def statsAndText(self, arr, mask=None):
        """Calculate some stats from an array and return them
        and some text.
        """
        numPoints = len(arr)
        if mask is not None:
            arr = arr[mask]
        med = np.nanmedian(arr)
        sigMad = nansigmaMad(arr)

        statsText = (
            "Median: {:0.2f}\n".format(med)
            + r"$\sigma_{MAD}$: "
            + "{:0.2f}\n".format(sigMad)
            + r"n$_{points}$: "
            + "{}".format(numPoints)
        )

        return med, sigMad, statsText

    def makePlot(
        self,
        data: KeyedData,
        plotInfo: Optional[Mapping[str, str]] = None,
        sumStats: Optional[Mapping] = None,
        **kwargs,
    ) -> Figure:
        """Make a skyPlot of the given data.

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
            and dec of the corners of the patch.

        Returns
        -------
        `pipeBase.Struct` containing:
            skyPlot : `matplotlib.figure.Figure`
                The resulting figure.

        Notes
        -----
        Expects the data to contain slightly different things
        depending on the types specified in plotTypes. This
        is handled automatically if you go through the pipetask
        framework but if you call this method separately then you
        need to make sure that data contains what the code is expecting.

        If stars is in the plot types given then it is expected that
        data contains: xStars, yStars, zStars and starStatMask.

        If galaxies is present: xGalaxies, yGalaxies, zGalaxies and
        galaxyStatsMask.

        If unknown is present: xUnknowns, yUnknowns, zUnknowns and
        unknownStatMask.

        If any is specified: x, y, z, statMask.

        These options are not exclusive and multiple can be specified
        and thus need to be present in data.

        Examples
        --------
        An example of the plot produced from this code is here:

        .. image:: /_static/analysis_tools/skyPlotExample.png

        For a detailed example of how to make a plot from the command line
        please see the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """

        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)

        if sumStats is None:
            sumStats = {}

        if plotInfo is None:
            plotInfo = {}

        # Make divergent colormaps for stars, galaxes and all the points
        blueGreen = mkColormap(["midnightblue", "lightcyan", "darkgreen"])
        redPurple = mkColormap(["indigo", "lemonchiffon", "firebrick"])
        orangeBlue = mkColormap(["darkOrange", "thistle", "midnightblue"])

        xCol = self.xAxisLabel
        yCol = self.yAxisLabel
        zCol = self.zAxisLabel  # noqa: F841

        toPlotList = []
        # For galaxies
        if "galaxies" in self.plotTypes:
            sortedArrs = sortAllArrays(
                [data["zGalaxies"], data["xGalaxies"], data["yGalaxies"], data["galaxyStatMask"]]
            )
            [colorValsGalaxies, xsGalaxies, ysGalaxies, statGalaxies] = sortedArrs
            statGalMed, statGalMad, galStatsText = self.statsAndText(colorValsGalaxies, mask=statGalaxies)
            # Add statistics
            bbox = dict(facecolor="lemonchiffon", alpha=0.5, edgecolor="none")
            # Check if plotting stars and galaxies, if so move the
            # text box so that both can be seen. Needs to be
            # > 2 becuase not being plotted points are assigned 0
            if len(self.plotTypes) > 2:
                boxLoc = 0.63
            else:
                boxLoc = 0.8
            ax.text(boxLoc, 0.91, galStatsText, transform=fig.transFigure, fontsize=8, bbox=bbox)
            toPlotList.append((xsGalaxies, ysGalaxies, colorValsGalaxies, redPurple, "Galaxies"))

        # For stars
        if "stars" in self.plotTypes:
            sortedArrs = sortAllArrays([data["zStars"], data["xStars"], data["yStars"], data["starStatMask"]])
            [colorValsStars, xsStars, ysStars, statStars] = sortedArrs
            statStarMed, statStarMad, starStatsText = self.statsAndText(colorValsStars, mask=statStars)
            # Add statistics
            bbox = dict(facecolor="paleturquoise", alpha=0.5, edgecolor="none")
            ax.text(0.8, 0.91, starStatsText, transform=fig.transFigure, fontsize=8, bbox=bbox)
            toPlotList.append((xsStars, ysStars, colorValsStars, blueGreen, "Stars"))

        # For unknowns
        if "unknown" in self.plotTypes:
            sortedArrs = sortAllArrays(
                [data["zUnknowns"], data["xUnknowns"], data["yUnknowns"], data["unknownStatMask"]]
            )
            [colorValsUnknowns, xsUnknowns, ysUnknowns, statUnknowns] = sortedArrs
            statUnknownMed, statUnknownMad, unknownStatsText = self.statsAndText(
                colorValsUnknowns, mask=statUnknowns
            )
            bbox = dict(facecolor="green", alpha=0.2, edgecolor="none")
            ax.text(0.8, 0.91, unknownStatsText, transform=fig.transFigure, fontsize=8, bbox=bbox)
            toPlotList.append((xsUnknowns, ysUnknowns, colorValsUnknowns, "viridis", "Unknown"))

        if "any" in self.plotTypes:
            sortedArrs = sortAllArrays([data["z"], data["x"], data["y"], data["statMask"]])
            [colorValsAny, xs, ys, statAny] = sortedArrs
            statAnyMed, statAnyMad, anyStatsText = self.statsAndText(colorValsAny, mask=statAny)
            bbox = dict(facecolor="purple", alpha=0.2, edgecolor="none")
            ax.text(0.8, 0.91, anyStatsText, transform=fig.transFigure, fontsize=8, bbox=bbox)
            toPlotList.append((xs, ys, colorValsAny, orangeBlue, "All"))

        # Corner plot of patches showing summary stat in each
        if self.plotOutlines:
            patches = []
            for dataId in sumStats.keys():
                (corners, _) = sumStats[dataId]
                ra = corners[0][0].asDegrees()
                dec = corners[0][1].asDegrees()
                xy = (ra, dec)
                width = corners[2][0].asDegrees() - ra
                height = corners[2][1].asDegrees() - dec
                patches.append(Rectangle(xy, width, height, alpha=0.3))
                ras = [ra.asDegrees() for (ra, dec) in corners]
                decs = [dec.asDegrees() for (ra, dec) in corners]
                ax.plot(ras + [ras[0]], decs + [decs[0]], "k", lw=0.5)
                cenX = ra + width / 2
                cenY = dec + height / 2
                if dataId == "tract":
                    minRa = np.min(ras)
                    minDec = np.min(decs)
                    maxRa = np.max(ras)
                    maxDec = np.max(decs)
                if dataId != "tract":
                    ax.annotate(
                        dataId,
                        (cenX, cenY),
                        color="k",
                        fontsize=5,
                        ha="center",
                        va="center",
                        path_effects=[pathEffects.withStroke(linewidth=2, foreground="w")],
                    )

        for i, (xs, ys, colorVals, cmap, label) in enumerate(toPlotList):
            finite = np.isfinite(xs) & np.isfinite(ys)
            xs = xs[finite]
            ys = ys[finite]
            n_xs = len(xs)
            # colorVal column is unusable so zero it out
            # This should be obvious on the plot
            if not any(np.isfinite(colorVals)):
                colorVals[:] = 0

            if n_xs < 5:
                continue
            if not self.plotOutlines or "tract" not in sumStats.keys():
                minRa = np.min(xs)
                maxRa = np.max(xs)
                minDec = np.min(ys)
                maxDec = np.max(ys)
                # Avoid identical end points which causes problems in binning
                if minRa == maxRa:
                    maxRa += 1e-5  # There is no reason to pick this number in particular
                if minDec == maxDec:
                    maxDec += 1e-5  # There is no reason to pick this number in particular

            plotOut = plotProjectionWithBinning(
                ax,
                xs,
                ys,
                colorVals,
                cmap,
                minRa,
                maxRa,
                minDec,
                maxDec,
                fixAroundZero=self.fixAroundZero,
                isSorted=True,
            )
            cax = fig.add_axes([0.87 + i * 0.04, 0.11, 0.04, 0.77])
            plt.colorbar(plotOut, cax=cax, extend="both")
            colorBarLabel = "{}: {}".format(self.zAxisLabel, label)
            text = cax.text(
                0.5,
                0.5,
                colorBarLabel,
                color="k",
                rotation="vertical",
                transform=cax.transAxes,
                ha="center",
                va="center",
                fontsize=10,
            )
            text.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])
            cax.tick_params(labelsize=7)

            if i == 0 and len(toPlotList) > 1:
                cax.yaxis.set_ticks_position("left")

        ax.set_xlabel(xCol)
        ax.set_ylabel(yCol)
        ax.tick_params(axis="x", labelrotation=25)
        ax.tick_params(labelsize=7)

        ax.set_aspect("equal")
        plt.draw()

        # Find some useful axis limits
        lenXs = [len(xs) for (xs, _, _, _, _) in toPlotList]
        if lenXs != [] and np.max(lenXs) > 1000:
            padRa = (maxRa - minRa) / 10
            padDec = (maxDec - minDec) / 10
            ax.set_xlim(maxRa + padRa, minRa - padRa)
            ax.set_ylim(minDec - padDec, maxDec + padDec)
        else:
            ax.invert_xaxis()

        # Add useful information to the plot
        plt.subplots_adjust(wspace=0.0, hspace=0.0, right=0.85)
        fig = plt.gcf()
        fig = addPlotInfo(fig, plotInfo)

        return fig
