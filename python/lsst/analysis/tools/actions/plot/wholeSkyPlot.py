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

__all__ = ("WholeSkyPlot",)

from collections.abc import Iterable
from typing import Mapping, Optional

import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import Config, Field, ListField
from matplotlib.collections import PatchCollection
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from ...math import nanSigmaMad
from .plotUtils import addPlotInfo, mkColormap


class WholeSkyPlot(PlotAction):
    """Plots the on sky distribution of a parameter.

    Plots the values of the parameter given for the z axis
    according to the positions given for x and y. Optimised
    for use with RA and Dec. Also calculates some basic
    statistics and includes those on the plot.

    The default axes limits and figure size were chosen to plot HSC PDR2.
    """

    keyBands = ListField[str](doc="Band for each metric, if any.")
    plotKeys = ListField[str](doc="Names of metrics to plot.")
    xAxisLabel = Field[str](doc="Label to use for the x axis.", default="RA (degrees)")
    yAxisLabel = Field[str](doc="Label to use for the y axis.", default="Dec (degrees)")
    figureSize = ListField[float](doc="Size of the figure.", default=[9.0, 3.5])
    colorBarRange = Field[float](
        doc="The multiplier for the color bar range. The max/min range values are: median +/- N * sigmaMad"
        ", where N is this config value.",
        default=3.0,
    )

    def getOutputNames(self, config: Config | None = None) -> Iterable[str]:
        return list(self.plotKeys)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        base.append(("corners", Vector))
        for key in self.plotKeys:
            base.append((key, Vector))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

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

    def _getAxesLimits(self, xs: list, ys: list) -> tuple(list, list):
        """Get the x and y axes limits in degrees.

        Parameters
        ----------
        xs : `list`
            X coordinates for the tracts to plot.
        ys : `list`
            Y coordinates for the tracts to plot.

        Returns
        -------
        xlim : `list`
            Minimun and maximum x axis values.
        ylim : `list`
            Minimun and maximum y axis values.
        """

        # Add some blank space on the edges of the plot.
        xlim = [np.nanmin(xs) - 5, np.nanmax(xs) + 5]
        ylim = [np.nanmin(ys) - 5, np.nanmax(ys) + 5]

        # Limit to only show real RA/Dec values.
        if xlim[0] < 0.0:
            xlim[0] = 0.0
        if xlim[1] > 360.0:
            xlim[1] = 360.0
        if ylim[0] < -90.0:
            ylim[0] = -90.0
        if ylim[1] > 90.0:
            ylim[1] = 90.0

        return (xlim, ylim)

    def _getMaxOutlierVals(self, multiplier: float, tracts: list, values: list, outlierInds: list) -> str:
        """Get the 5 largest outlier values in a string.

        Parameters
        ----------
        multiplier : `float`
            Select values whose absolute value is > multiplier * sigmaMAD.
        tracts : `list`
            All the tracts.
        values : `list`
            All the metric values.
        outlierInds : `list`
            Indicies of outlier values.

        Returns
        -------
        text : `str`
            A string containing the 5 tracts with the largest outlier values.
        """
        text = f"Tracts with |value| > {multiplier}" + r"$\sigma_{MAD}$" + ": "
        if len(outlierInds) > 0:
            outlierValues = np.array(values)[outlierInds]
            outlierTracts = np.array(tracts)[outlierInds]
            # Sort values in descending (-) absolute value order discounting
            # NaNs.
            maxInds = np.argsort(-np.abs(outlierValues))
            # Show up to five values on the plot.
            for ind in maxInds[:5]:
                val = outlierValues[ind]
                tract = outlierTracts[ind]
                text += f"{tract}, {val:.3}; "
            # Remove the final trailing comma and whitespace.
            text = text[:-2]
        else:
            text += "None"

        return text

    def makePlot(
        self,
        data: KeyedData,
        plotInfo: Optional[Mapping[str, str]] = None,
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

        Returns
        -------
        results : `dict`
            A dictionary containing the resulting figures.


        Examples
        --------
        An example of the plot produced from this code is here:

        .. image:: /_static/analysis_tools/wholeSkyPlotExample.png

        For a detailed example of how to make a plot from the command line
        please see the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """
        # Make a divergent colormap.
        blueGreen = mkColormap(["midnightblue", "lightcyan", "darkgreen"])

        results = {}
        for key, band in zip(self.plotKeys, self.keyBands):

            if plotInfo is None:
                plotInfo = {}

            # Create patches using the corners of each tract.
            patches = []
            colBarVals = []
            tracts = []
            ras = []
            decs = []
            for i, corners in enumerate(data["corners"]):
                patches.append(Polygon(corners, closed=True))
                colBarVals.append(data[key][i])
                tracts.append(data["tract"][i])
                ras.append(corners[0][0])
                decs.append(corners[0][1])

            # Setup figure.
            fig, ax = plt.subplots(1, 1, figsize=self.figureSize, dpi=500)
            xlim, ylim = self._getAxesLimits(ras, decs)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(self.xAxisLabel)
            ax.set_ylabel(self.yAxisLabel)
            ax.invert_xaxis()

            # Add colored patches showing tract metric values.
            patchCollection = PatchCollection(patches, cmap=blueGreen)
            ax.add_collection(patchCollection)

            # Define color bar range.
            med = np.nanmedian(colBarVals)
            sigmaMad = nanSigmaMad(colBarVals)
            vmin = med - self.colorBarRange * sigmaMad
            vmax = med + self.colorBarRange * sigmaMad

            # Note tracts with metrics outside (vmin, vmax) as outliers.
            outlierInds = np.where((colBarVals < vmin) | (colBarVals > vmax))[0]

            # Initialize legend handles.
            handles = []

            # Plot the outlier patches.
            outlierPatches = []
            if len(outlierInds) > 0:
                for ind in outlierInds:
                    outlierPatches.append(patches[ind])
                outlierPatchCollection = PatchCollection(
                    outlierPatches,
                    cmap=blueGreen,
                    facecolors="none",
                    edgecolors="r",
                    linewidths=0.5,
                    zorder=100,
                )
                ax.add_collection(outlierPatchCollection)
                # Add legend information.
                outlierPatch = Patch(
                    facecolor="none",
                    edgecolor="r",
                    linewidth=0.5,
                    label="Outlier",
                )
                handles.append(outlierPatch)

            # Plot tracts with NaN metric values.
            nanInds = np.where(~np.isfinite(colBarVals))[0]
            nanPatches = []
            if len(nanInds) > 0:
                for ind in nanInds:
                    nanPatches.append(patches[ind])
                nanPatchCollection = PatchCollection(
                    nanPatches,
                    cmap=None,
                    facecolors="white",
                    edgecolors="grey",
                    linestyles="dotted",
                    linewidths=0.5,
                    zorder=10,
                )
                ax.add_collection(nanPatchCollection)
                # Add legend information.
                nanPatch = Patch(
                    facecolor="white",
                    edgecolor="grey",
                    linestyle="dotted",
                    linewidth=0.5,
                    label="NaN",
                )
                handles.append(nanPatch)

            if len(handles) > 0:
                fig.legend(handles=handles)

            # Add text boxes to show the number of tracts, number of NaNs,
            # median, sigma MAD, and the five largest outlier values.
            outlierText = self._getMaxOutlierVals(self.colorBarRange, tracts, colBarVals, outlierInds)
            multiplier = 3.5 / self.figureSize[1]
            verticalSpacing = 0.028 * multiplier
            fig.text(
                0.01,
                0.01 + 3 * verticalSpacing,
                f"Num tracts: {len(tracts)}",
                transform=fig.transFigure,
                fontsize=8,
                alpha=0.7,
            )
            fig.text(
                0.01,
                0.01 + 2 * verticalSpacing,
                f"Num nans: {len(nanInds)}",
                transform=fig.transFigure,
                fontsize=8,
                alpha=0.7,
            )
            fig.text(
                0.01,
                0.01 + verticalSpacing,
                f"Median: {med:.3f}; " + r"$\sigma_{MAD}$" + f": {sigmaMad:.3f}",
                transform=fig.transFigure,
                fontsize=8,
                alpha=0.7,
            )
            fig.text(0.01, 0.01, outlierText, transform=fig.transFigure, fontsize=8, alpha=0.7)

            # Truncate the color range to (vmin, vmax).
            colorBarVals = np.clip(np.array(colBarVals), vmin, vmax)
            patchCollection.set_array(colorBarVals)
            # Make the color bar with a metric label.
            cbar = plt.colorbar(
                patchCollection,
                ax=ax,
                shrink=0.7,
                extend="both",
                location="top",
                orientation="horizontal",
            )
            cbarText = key
            if data[key].unit is not None:
                unit = data[key].unit.to_string()
                cbarText += f" ({unit})"
            text = cbar.ax.text(
                0.5,
                0.5,
                cbarText,
                transform=cbar.ax.transAxes,
                ha="center",
                va="center",
                fontsize=10,
                zorder=100,
            )
            text.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])

            # Finalize plot appearance.
            ax.grid()
            ax.set_axisbelow(True)
            ax.set_aspect("equal")
            fig = plt.gcf()
            plotInfo["bands"] = band
            fig = addPlotInfo(fig, plotInfo)
            plt.subplots_adjust(left=0.08, right=0.97, top=0.8, bottom=0.17, wspace=0.35)

            results[key] = fig
        return results
