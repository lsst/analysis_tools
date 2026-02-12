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

import importlib.resources as importResources
import json
from typing import Mapping, Optional

import matplotlib.patheffects as pathEffects
import numpy as np
import yaml
from matplotlib import gridspec
from matplotlib.collections import PatchCollection
from matplotlib.colors import CenteredNorm
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon

import lsst.analysis.tools
from lsst.pex.config import ChoiceField, Field, ListField
from lsst.utils.plotting import (
    accent_color,
    divergent_cmap,
    make_figure,
    set_rubin_plotstyle,
    stars_cmap,
    stars_color,
)

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from ...math import nanSigmaMad
from ...utils import getTractCorners
from .plotUtils import addPlotInfo


class AnnotatedFigure(Figure):
    metadata: dict


class WholeSkyPlot(PlotAction):
    """Plots the on sky distribution of a parameter.

    Plots the values of the parameter given for the z axis
    according to the positions given for x and y. Optimised
    for use with RA and Dec. Also calculates some basic
    statistics and includes those on the plot.

    The default axes limits and figure size were chosen to plot HSC PDR2.
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis.", default="RA (degrees)")
    yAxisLabel = Field[str](doc="Label to use for the y axis.", default="Dec (degrees)")
    zAxisLabel = Field[str](doc="Label to use for the z axis.", default="")
    autoAxesLimits = Field[bool](doc="Find axes limits automatically.", default=True)
    xLimits = ListField[float](doc="Plotting limits for the x axis.", default=[-5.0, 365.0])
    yLimits = ListField[float](doc="Plotting limits for the y axis.", default=[-10.0, 60.0])
    autoAxesLimits = Field[bool](doc="Find axes limits automatically.", default=True)
    colorBarMin = Field[float](doc="The minimum value of the color bar.", optional=True)
    colorBarMax = Field[float](doc="The minimum value of the color bar.", optional=True)
    colorBarRange = Field[float](
        doc="The multiplier for the color bar range. The max/min range values are: median +/- N * sigmaMad"
        ", where N is this config value.",
        default=3.0,
    )
    colorMapType = ChoiceField[str](
        doc="Type of color map to use for the color bar. Options: sequential, divergent, userDefined.",
        allowed={cmType: cmType for cmType in ("sequential", "divergent")},
        default="divergent",
    )
    colorMap = ListField[str](
        doc="List of hexidecimal colors for a user-defined color map.",
        optional=True,
    )
    showOutliers = Field[bool](
        doc="Show the outliers on the plot. "
        "Outliers are values whose absolute value is > colorBarRange * sigmaMAD.",
        default=True,
    )
    showNaNs = Field[bool](doc="Show the NaNs on the plot.", default=True)
    labelTracts = Field[bool](doc="Label the tracts.", default=False)

    addThresholds = Field[bool](
        doc="Read in the predefined thresholds and indicate them on the histogram.",
        default=True,
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        base.append(("z", Vector))
        base.append(("tract", Vector))
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
            A string containing the 10 tracts with the largest outlier values.
        """
        if self.addThresholds:
            text = "Tracts with value outside thresholds: "
        else:
            text = f"Tracts with |value| > {multiplier}" + r"$\sigma_{MAD}$" + ": "
        if len(outlierInds) > 0:
            outlierValues = np.array(values)[outlierInds]
            outlierTracts = np.array(tracts)[outlierInds]
            # Sort values in descending (-) absolute value order discounting
            # NaNs.
            maxInds = np.argsort(-np.abs(outlierValues))
            # Show up to ten values on the plot.
            for ind in maxInds[:10]:
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
    ) -> AnnotatedFigure:
        """Make a WholeSkyPlot of the given data.

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
        `pipeBase.Struct` containing:
            skyPlot : `matplotlib.figure.Figure`
                The resulting figure.


        Examples
        --------
        An example of the plot produced from this code is here:

        .. image:: /_static/analysis_tools/wholeSkyPlotExample.png

        For a detailed example of how to make a plot from the command line
        please see the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """
        skymap = kwargs["skymap"]
        if plotInfo is None:
            plotInfo = {}

        if self.addThresholds:
            metricThresholdFile = importResources.read_text(lsst.analysis.tools, "metricInformation.yaml")
            metricDefs = yaml.safe_load(metricThresholdFile)

        # Prevent Bands in the plot info showing a list of bands.
        # If bands is a list, it implies that parameterizedBand=False,
        # and that the metric is not band-specific.
        if "bands" in plotInfo:
            if isinstance(plotInfo["bands"], list):
                plotInfo["bands"] = "N/A"

        colorMap = self.colorMap
        match self.colorMapType:
            case "sequential":
                if colorMap is None:
                    colorMap = stars_cmap()
                    outlierColor = "red"
                norm = None
            case "divergent":
                if colorMap is None:
                    colorMap = divergent_cmap()
                    outlierColor = "fuchsia"
                norm = CenteredNorm()

        # Create patches using the corners of each tract.
        patches = []
        colBarVals = []
        tracts = []
        ras = []
        decs = []
        mid_ras = []
        mid_decs = []
        for i, tract in enumerate(data["tract"]):
            corners = getTractCorners(skymap, tract)
            patches.append(Polygon(corners, closed=True))
            colBarVals.append(data["z"][i])
            tracts.append(tract)
            ras.append(corners[0][0])
            decs.append(corners[0][1])
            mid_ras.append((corners[0][0] + corners[1][0]) / 2)
            mid_decs.append((corners[0][1] + corners[2][1]) / 2)

        # Setup figure.
        fig: AnnotatedFigure = make_figure(dpi=300, figsize=(12, 3.5))
        set_rubin_plotstyle()
        gs = gridspec.GridSpec(1, 4)
        ax = fig.add_subplot(gs[:3])
        # Add colored patches showing tract metric values.
        patchCollection = PatchCollection(patches, cmap=colorMap, norm=norm)
        ax.add_collection(patchCollection)

        # Define color bar range.
        if np.sum(np.isfinite(colBarVals)) > 0:
            med = np.nanmedian(colBarVals)
        else:
            med = np.nan
        sigmaMad = nanSigmaMad(colBarVals)
        if self.colorBarMin is not None:
            vmin = np.float64(self.colorBarMin)
        else:
            vmin = med - self.colorBarRange * sigmaMad
        if self.colorBarMax is not None:
            vmax = np.float64(self.colorBarMax)
        else:
            vmax = med + self.colorBarRange * sigmaMad

        dataName = self.zAxisLabel.format_map(kwargs)
        colBarVals = np.array(colBarVals)
        if self.addThresholds and dataName in metricDefs:
            if "lowThreshold" in metricDefs[dataName].keys():
                lowThreshold = metricDefs[dataName]["lowThreshold"]
            else:
                lowThreshold = np.nan
            if "highThreshold" in metricDefs[dataName].keys():
                highThreshold = metricDefs[dataName]["highThreshold"]
            else:
                highThreshold = np.nan
            outlierInds = np.where((colBarVals < lowThreshold) | (colBarVals > highThreshold))[0]
        else:
            # Note tracts with metrics outside (vmin, vmax) as outliers.
            outlierInds = np.where((colBarVals < vmin) | (colBarVals > vmax))[0]

        # Initialize legend handles.
        handles = []

        if self.showOutliers:
            # Plot the outlier patches.
            outlierPatches = []
            if len(outlierInds) > 0:
                for ind in outlierInds:
                    outlierPatches.append(patches[ind])
                outlierPatchCollection = PatchCollection(
                    outlierPatches,
                    cmap=colorMap,
                    norm=norm,
                    facecolors="none",
                    edgecolors=outlierColor,
                    linewidths=0.5,
                    zorder=100,
                )
                ax.add_collection(outlierPatchCollection)
                # Add legend information.
                outlierPatch = Patch(
                    facecolor="none",
                    edgecolor=outlierColor,
                    linewidth=0.5,
                    label="Outlier",
                )
                handles.append(outlierPatch)

        if self.showNaNs:
            # Plot tracts with NaN metric values.
            nanInds = np.where(~np.isfinite(colBarVals))[0]
            nanPatches = []
            if len(nanInds) > 0:
                for ind in nanInds:
                    nanPatches.append(patches[ind])
                nanPatchCollection = PatchCollection(
                    nanPatches,
                    cmap=None,
                    norm=norm,
                    facecolors="white",
                    edgecolors="grey",
                    linestyles="dotted",
                    linewidths=0.5,
                    zorder=100,
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

        if self.labelTracts:
            # Label the tracts
            for i, tract in enumerate(tracts):
                ax.text(
                    mid_ras[i],
                    mid_decs[i],
                    f"{tract}",
                    ha="center",
                    va="center",
                    fontsize=2,
                    alpha=0.7,
                    zorder=100,
                )

        ax.set_aspect("equal")
        axPos = ax.get_position()
        ax1 = fig.add_axes([0.73, 0.25, 0.20, 0.47])

        if np.sum(np.isfinite(data["z"])) > 0:
            ax1.hist(data["z"], bins=len(data["z"] / 10), color=stars_color(), histtype="step")
        else:
            ax1.text(0.5, 0.5, "Data all NaN/Inf")
        ax1.set_xlabel("Metric Values")
        ax1.set_ylabel("Number")
        ax1.yaxis.set_label_position("right")
        ax1.yaxis.tick_right()

        if self.addThresholds and dataName in metricDefs:
            # Check the thresholds are finite and set them to
            # the min/max of the data if they aren't to calculate
            # the x range of the plot
            if np.isfinite(lowThreshold):
                ax1.axvline(lowThreshold, color=accent_color())
            else:
                lowThreshold = np.nanmin(colBarVals)
            if np.isfinite(highThreshold):
                ax1.axvline(highThreshold, color=accent_color())
            else:
                highThreshold = np.nanmax(colBarVals)

            widthThreshold = highThreshold - lowThreshold
            upperLim = highThreshold + 0.5 * widthThreshold
            lowerLim = lowThreshold - 0.5 * widthThreshold
            ax1.set_xlim(lowerLim, upperLim)
            numOutside = np.sum(((data["z"] > upperLim) | (data["z"] < lowerLim)))
            ax1.set_title("Outside plot limits: " + str(numOutside))

        else:
            if vmin != vmax and np.isfinite(vmin) and np.isfinite(vmax):
                ax1.set_xlim(vmin, vmax)

        if self.autoAxesLimits:
            xlim, ylim = self._getAxesLimits(ras, decs)
        else:
            xlim, ylim = self.xLimits, self.yLimits
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)
        ax.invert_xaxis()

        if self.showOutliers:
            # Add text boxes to show the number of tracts, number of NaNs,
            # median, sigma MAD, and the five largest outlier values.
            outlierText = self._getMaxOutlierVals(self.colorBarRange, tracts, colBarVals, outlierInds)
        # Make vertical text spacing readable for different figure sizes.
        multiplier = 3.5 / fig.get_size_inches()[1]
        verticalSpacing = 0.028 * multiplier
        fig.text(
            0.01,
            0.01 + 3 * verticalSpacing,
            f"Num tracts: {len(tracts)}",
            transform=fig.transFigure,
            fontsize=8,
            alpha=0.7,
        )
        if self.showNaNs:
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
        if self.showOutliers:
            fig.text(0.01, 0.01, outlierText, transform=fig.transFigure, fontsize=8, alpha=0.7)

        # Truncate the color range to (vmin, vmax).
        if vmin != vmax and np.isfinite(vmin) and np.isfinite(vmax):
            colBarVals = np.clip(np.array(colBarVals), vmin, vmax)
        patchCollection.set_array(colBarVals)
        # Make the color bar with a metric label.
        axPos = ax.get_position()
        cax = fig.add_axes([0.084, axPos.y1 + 0.02, 0.62, 0.07])
        fig.colorbar(
            patchCollection,
            cax=cax,
            shrink=0.7,
            extend="both",
            location="top",
            orientation="horizontal",
        )
        cbarText = "Metric Values"

        text = cax.text(
            0.5,
            0.5,
            cbarText,
            transform=cax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
            zorder=100,
        )
        text.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])

        # Finalize plot appearance.
        ax.grid()
        ax.set_axisbelow(True)
        addPlotInfo(fig, plotInfo)
        fig.subplots_adjust(left=0.08, right=0.92, top=0.8, bottom=0.17, wspace=0.05)
        titleText = self.zAxisLabel.format_map(kwargs)
        if "zUnit" in data and data["zUnit"] != "":
            titleText += f" ({data['zUnit']})"
        fig.suptitle("Metric: " + titleText, fontsize=20)

        # This saves metadata in the PNG that allows the plot-navigator
        # to provide tract numbers and metric values on mouseover.
        #
        # PNG metadata is a set of string keys and string values.
        # The WholeSkyPlot stores two keys:
        # - label: the string describing the regions ('tract')
        # - boxes, JSON string of a list of per-region dictionaries,
        # where each dictionary has fields:
        # - min_x, max_x, min_y, max_y for the pixel coordinates of
        #   the four corners of the region
        # - id: the identifier of the region (e.g. tract number)
        # - value: the region's metric, as a string.
        #
        def make_patch_md(patch, id_field, value, ax):
            path = ax.transData.transform_path(patch.get_path())
            x_path = [int(x) for x in path.vertices[:, 0].tolist()]
            y_path = [int(y) for y in path.vertices[:, 1].tolist()]
            return {
                "min_x": min(x_path),
                "max_x": max(x_path),
                "min_y": min(y_path),
                "max_y": max(y_path),
                "id": f"{id_field}",
                "value": f"{value:.3}",
            }

        # After ax.set_aspect(), the figure needs to be drawn for the axes
        # transformations to be updated to the right values.
        fig.canvas.draw_idle()

        patch_coordinate_entries = [
            make_patch_md(patch, tract, value, ax)
            for (patch, tract, value) in zip(patches, tracts, colBarVals)
        ]

        fig.metadata = {"label": "Tract", "boxes": json.dumps(patch_coordinate_entries)}

        return fig
