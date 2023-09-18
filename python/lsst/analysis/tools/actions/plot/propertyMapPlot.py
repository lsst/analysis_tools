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

__all__ = ("PropertyMapPlot",)

import logging
from typing import Iterable, Mapping, Union

import lsst.pex.config as pexConfig
import matplotlib.patheffects as mpl_path_effects
import matplotlib.pyplot as plt
import numpy as np
import skyproj
from healsparse.healSparseMap import HealSparseMap
from lsst.analysis.tools.tasks.propertyMapTractAnalysis import PropertyMapTractAnalysisConfig
from lsst.skymap.tractInfo import ExplicitTractInfo
from matplotlib.figure import Figure
from matplotlib.legend_handler import HandlerTuple

from ...interfaces import KeyedData, PlotAction

_LOG = logging.getLogger(__name__)


class CustomHandler(HandlerTuple):
    """Custom legend handler to overlay multiple patches for a single
    legend entry.

    This handler class inherits from `HandlerTuple` and is designed to
    handle cases where multiple artists (e.g., patches) need to be overlaid
    on top of each other in a single legend entry, as opposed to
    side-by-side which is the default behavior of `HandlerTuple`.

    Methods
    -------
    create_artists:
        Override the `create_artists` method of `HandlerTuple` to modify
        the positioning of the artists so that they overlay directly on top
        of one another in the legend.

    Example
    -------
    # Plot some data.
    line, = ax.plot(x, y, label='Sample Line')

    # Use CustomHandler for overlaid patches and also include the regular
    # line legend if desired.
    handles = [(patch1, patch2), line]
    labels = ['Overlaid Patches', line.get_label()]
    leg = ax.legend(
        handles, labels, handler_map={tuple: CustomHandler()}, loc="best"
    )
    """

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        artists = HandlerTuple.create_artists(
            self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
        )
        # Overlay the two patches.
        for a in artists:
            a.set_transform(trans)
        return artists


class PropertyMapPlot(PlotAction):
    plotName = pexConfig.Field[str](doc="The name for the plotting task.", optional=True)

    def __call__(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        plotConfig: PropertyMapTractAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
        **kwargs,
    ) -> Mapping[str, Figure]:
        self._validateInput(data, tractInfo, plotConfig, plotInfo)
        return self.makePlot(data, tractInfo, plotConfig, plotInfo)

    def _validateInput(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        plotConfig: PropertyMapTractAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
    ) -> None:
        """Validate the input data."""

        if not isinstance(tractInfo, ExplicitTractInfo):
            raise TypeError(f"Input `tractInfo` type must be {ExplicitTractInfo} not {type(tractInfo)}.")

        if not isinstance(plotConfig, PropertyMapTractAnalysisConfig):
            raise TypeError(
                "`plotConfig` must be a "
                "`lsst.analysis.tools.tasks.propertyMapTractAnalysis.PropertyMapTractAnalysisConfig`."
            )

        if not isinstance(plotInfo, dict):
            raise TypeError("`plotConfig` must be a dictionary.")

        zoomFactors = plotConfig.zoomFactors
        isListOfFloats = isinstance(zoomFactors, pexConfig.listField.List) and all(
            isinstance(zf, float) for zf in zoomFactors
        )
        if not (isListOfFloats and len(zoomFactors) == 2) or any(zf <= 1 for zf in zoomFactors):
            raise TypeError(
                "`zoomFactors` must be a two-element `lsst.pex.config.listField.List` of floats > 1."
            )

        for prop, propConfig in plotConfig.properties.items():
            if not isinstance(propConfig.nBinsHist, int) or propConfig.nBinsHist <= 0:
                raise ValueError(
                    f"`nBinsHist` for property `{prop}` must be a positive integer, not "
                    f"{propConfig.nBinsHist}."
                )

        # Identify any invalid entries in `data`.
        invalidEntries = {
            key: pytype
            for key, pytype in {k: v.ref.datasetType.storageClass.pytype for k, v in data.items()}.items()
            if pytype != HealSparseMap
        }

        # If any invalid entries are found, raise a TypeError with details.
        if invalidEntries:
            errorMessage = "; ".join(
                f"`{key}` should be {HealSparseMap}, got {type_}" for key, type_ in invalidEntries.items()
            )
            raise TypeError(f"Invalid input types found in `data`: {errorMessage}")

    def addPlotInfo(
        self,
        fig: Figure,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
        mapName: Mapping[str, str],
    ) -> Figure:
        """Add useful information to the plot.

        Parameters
        ----------
        fig : `matplotlib.figure.Figure`
            The figure to add the information to.
        plotInfo : `dict`
            A dictionary of the plot information.
        mapName : `str`
            The name of the map being plotted.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The figure with the information added.
        """

        run = plotInfo["run"]
        tableType = f"\nTable: {plotInfo['tableNames'][mapName]}"

        dataIdText = f"Tract: {plotInfo['tract']}, Band: {plotInfo['band']}"
        mapText = (
            f", Property: {plotInfo['property']}, "
            f"Operation: {plotInfo['operation']}, "
            f"Coadd: {plotInfo['coaddName']}"
        )
        geomText = f", Valid area: {plotInfo['valid_area']:.2f} sq. deg., " f"NSIDE: {plotInfo['nside']}"
        infoText = f"\n{dataIdText}{mapText}"

        fig.text(
            0.04,
            0.965,
            f'{plotInfo["plotName"]}: {plotInfo["property"].title().replace("Psf", "PSF")}',
            fontsize=19,
            transform=fig.transFigure,
            ha="left",
            va="top",
        )
        t = fig.text(
            0.04,
            0.942,
            f"{run}{tableType}{geomText}{infoText}",
            fontsize=15,
            transform=fig.transFigure,
            alpha=0.6,
            ha="left",
            va="top",
        )
        t.set_linespacing(1.4)

        return fig

    def makePlot(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        plotConfig: PropertyMapTractAnalysisConfig,
        plotInfo: Mapping[str, Union[Mapping[str, str], str, int]],
    ) -> Mapping[str, Figure]:
        """Make the survey property map plot.

        Parameters
        ----------
        data : `KeyedData`
            The HealSparseMap to plot the points from.
        tractInfo: `lsst.skymap.tractInfo.ExplicitTractInfo`
            The tract info object.
        plotConfig :
            `lsst.analysis.tools.tasks.propertyMapTractAnalysis.
            PropertyMapTractAnalysisConfig`
            The configuration for the plot.
        plotInfo : `dict`
            A dictionary of information about the data being plotted.

        Returns
        -------
        figDict : `dict` [`~matplotlib.figure.Figure`]
            The resulting figures.
        """

        figDict: dict[str, Figure] = {}

        # Plotting customization.
        colorbarTickLabelSize = 14
        colorBarAspect = 16
        rcparams = {
            "axes.labelsize": 18,
            "axes.linewidth": 1.8,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
        }
        zoomFactors = plotConfig.zoomFactors

        # Muted green for the full map, and muted red and blue for the two
        # zoomed-in maps.
        histColors = ["#265D40", "#8B0000", "#00008B"]

        with plt.rc_context(rcparams):
            for mapName, deferredDatasetHandle in data.items():
                mapData = deferredDatasetHandle.get()
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 16))

                # Reduce whitespace but leave some room at the top for
                # `plotInfo`.
                plt.subplots_adjust(left=0.064, right=0.96, top=0.855, bottom=0.07, wspace=0.18, hspace=0.24)

                # Get the values for the valid pixels of the full tract.
                values = mapData[mapData.valid_pixels]
                goodValues = np.isfinite(values)
                values = values[goodValues]  # As a precaution.

                # Make a concise human-readable label for the plot.
                plotInfo["coaddName"] = mapName.split("Coadd_")[0]
                plotInfo["operation"] = self.getLongestSuffixMatch(
                    mapName, ["min", "max", "mean", "weighted_mean", "sum"]
                ).replace("_", " ")
                propertyName = mapName[
                    len(f"{plotInfo['coaddName']}Coadd_") : -len(plotInfo["operation"])
                ].strip("_")
                plotInfo["property"] = propertyName.replace("_", " ")

                nBinsHist = plotConfig.properties[propertyName].nBinsHist
                fullExtent = None
                zoomIdx = []
                for ax, zoom, zoomFactor, histColor in zip(
                    [ax1, ax3, ax4], [True, False, False], [None, *zoomFactors], histColors
                ):
                    extent = self.getZoomedExtent(fullExtent, zoomFactor)
                    sp = skyproj.GnomonicSkyproj(
                        ax=ax,
                        lon_0=tractInfo.ctr_coord.getRa().asDegrees(),
                        lat_0=tractInfo.ctr_coord.getDec().asDegrees(),
                        extent=extent,
                        rcparams=rcparams,
                    )
                    sp.draw_hspmap(mapData, zoom=zoom)
                    sp.set_xlabel("RA")
                    sp.set_ylabel("Dec")
                    cbar = sp.draw_colorbar(location="right", fraction=0.15, aspect=colorBarAspect, pad=0)
                    cbar.ax.tick_params(labelsize=colorbarTickLabelSize)
                    cbarText = (
                        "Full Tract" if zoomFactor is None else f"{self.prettyPrintFloat(zoomFactor)}x Zoom"
                    )
                    self.addTextToColorbar(cbar, cbarText, color=histColor)
                    if zoomFactor is None:
                        # Save the skyproj object of the full-tract plot.
                        # Will be used in drawing zoom rectangles etc.
                        spf = sp
                        # Get the extent of the full tract.
                        fullExtent = spf.get_extent()
                    else:
                        # Create a rectangle for the zoomed-in region.
                        x0, x1, y0, y1 = extent
                        for c, ls, lw in zip(["white", histColor], ["solid", "dashed"], [3.5, 1.5]):
                            spf.draw_polygon(
                                [x0, x0, x1, x1],
                                [y0, y1, y1, y0],
                                facecolor="none",
                                edgecolor=c,
                                linestyle=ls,
                                linewidth=lw,
                                alpha=0.8,
                            )
                        zoomText = spf.ax.text(
                            (x0 + x1) / 2,
                            y0,
                            f"{self.prettyPrintFloat(zoomFactor)}x",
                            color=histColor,
                            fontsize=14,
                            fontweight="bold",
                            alpha=0.8,
                            ha="center",
                            va="bottom",
                        )
                        # Add a distinct outline around the text for better
                        # visibility in various backgrounds.
                        zoomText.set_path_effects(
                            [
                                mpl_path_effects.Stroke(linewidth=4, foreground="white", alpha=0.8),
                                mpl_path_effects.Normal(),
                            ]
                        )
                        # Get the indices of pixels in the zoomed-in region.
                        pos = mapData.valid_pixels_pos()
                        # Reversed axes consideration.
                        xmin, xmax = sorted([x0, x1])
                        idx = (pos[0] > xmin) & (pos[0] < xmax) & (pos[1] > y0) & (pos[1] < y1)
                        zoomIdx.append(idx[goodValues])

                # Calculate weights for each bin to ensure that the peak of the
                # histogram reaches 1.
                weights = np.ones_like(values) / np.histogram(values, bins=nBinsHist)[0].max()

                # Compute full-tract histogram and get its bins.
                # NOTE: `exposure_time` histograms are quantized and look more
                # bar-like, so they look better with fewer bins.
                bins = ax2.hist(
                    values,
                    bins=nBinsHist,
                    label="Full Tract",
                    color=histColors[0],
                    weights=weights,
                    alpha=0.7,
                )[1]

                # Align the histogram (top right panel) with the skyproj plots.
                pos1 = spf.ax.get_position()  # Top left.
                pos4 = sp.ax.get_position()  # Bottom right.
                cbarWidth = cbar.ax.get_position().height / colorBarAspect
                # NOTE: cbarWidth != cbar.ax.get_position().width
                ax2.set_position([pos4.x0, pos1.y0, pos4.width + cbarWidth, pos1.height])

                # Overplot the histograms for the zoomed-in plots.
                for zoomFactor, zidx, color, linestyle, hatch in zip(
                    zoomFactors, zoomIdx, histColors[1:], ["solid", "dotted"], ["//", "xxxx"]
                ):
                    weights = np.ones_like(values[zidx]) / np.histogram(values[zidx], bins=bins)[0].max()
                    histLabel = f"{self.prettyPrintFloat(zoomFactor)}x Zoom"
                    histValues = ax2.hist(
                        values[zidx],
                        bins=bins,
                        label=histLabel,
                        color=color,
                        weights=weights,
                        histtype="step",
                        linewidth=2,
                        linestyle=linestyle,
                        alpha=0.6,
                    )[0]
                    # Fill the area under the step.
                    ax2.fill_between(
                        bins[:-1],
                        histValues,
                        step="post",
                        color=color,
                        alpha=0.2,
                        hatch=hatch,
                        label="hidden",
                    )

                # Set labels and legend.
                ax2.set_xlabel(plotInfo["property"].title().replace("Psf", "PSF"))
                ax2.set_ylabel("Normalized Count")

                # Get handles and labels from the axes.
                handles, labels = ax2.get_legend_handles_labels()

                # Add a legend with custom handler that combines the handle
                # pairs for the zoomed-in cases.
                handles = [handles[0], (handles[1], handles[2]), (handles[3], handles[4])]
                while "hidden" in labels:
                    labels.remove("hidden")
                legend = ax2.legend(
                    handles,
                    labels,
                    handler_map={tuple: CustomHandler()},
                    loc="best",
                    frameon=False,
                    fontsize=15,
                )

                for line, text in zip(handles, legend.get_texts()):
                    if isinstance(line, tuple):
                        # Use the first handle to get the color.
                        line = line[0]
                    color = line.get_edgecolor() if line.get_facecolor()[-1] == 0 else line.get_facecolor()
                    text.set_color(color)

                # Add extra info to `plotInfo`.
                plotInfo["nside"] = mapData.nside_sparse
                plotInfo["valid_area"] = mapData.get_valid_area()

                # Add useful information to the plot.
                figDict[mapName] = self.addPlotInfo(fig, plotInfo, mapName)

                _LOG.info(
                    f"Made property map plot for dataset type {mapName}, tract: {plotInfo['tract']}, "
                    f"band: '{plotInfo['band']}'."
                )

        return figDict

    def getOutputNames(self, config=None) -> Iterable[str]:
        # Docstring inherited.

        # Names needed for making corresponding output connections for the maps
        # that are configured for this task.
        outputNames: tuple[str] = ()
        for propertyName in config.properties:
            coaddName = config.properties[propertyName].coaddName
            for operationName in config.properties[propertyName].operations:
                outputNames += (f"{coaddName}Coadd_{propertyName}_{operationName}",)

        return outputNames

    @staticmethod
    def getZoomedExtent(fullExtent, n):
        """Get zoomed extent centered on the original full plot.

        Parameters
        ----------
        fullExtent : `tuple` of `float`
            The full extent defined by (lon_min, lon_max, lat_min, lat_max):

            * ``lon_min``
                Minimum longitude of the original extent (`float`).
            * ``"lon_max"``
                Maximum longitude of the original extent (`float`).
            * ``lat_min``
                Minimum latitude of the original extent (`float`).
            * ``"lat_max"``
                Maximum latitude of the original extent (`float`).

        n : `float`, optional
            Zoom factor; for instance, n=2 means zooming in 2 times at the
            center. If None, the function returns None.

        Returns
        -------
        Results : `tuple` [`float`]
            New extent as (new_lon_min, new_lon_max, new_lat_min, new_lat_max).
        """
        if n is None:
            return None
        lon_min, lon_max, lat_min, lat_max = fullExtent
        lon_center, lat_center = (lon_min + lon_max) / 2, (lat_min + lat_max) / 2
        half_lon = (lon_max - lon_min) * np.cos(np.radians(lat_center)) / (2 * n)
        half_lat = (lat_max - lat_min) / (2 * n)
        return lon_center - half_lon, lon_center + half_lon, lat_center - half_lat, lat_center + half_lat

    @staticmethod
    def prettyPrintFloat(n):
        if n.is_integer():
            return str(int(n))
        return str(n)

    @staticmethod
    def addTextToColorbar(
        cb, text, orientation="vertical", color="black", fontsize=14, fontweight="bold", alpha=0.8
    ):
        """Helper method to add text inside the horizontal colorbar.

        Parameters
        ----------
        cb : `matplotlib.colorbar.Colorbar`
            The colorbar object.
        text : `str`
            The text to add.
        orientation : `str`, optional
            The orientation of the colorbar. Can be either "vertical" or
            "horizontal".
        fontsize : `int`, optional
            The fontsize of the text.
        fontweight : `str`, optional
            The fontweight of the text.
        alpha : `float`, optional
            The alpha value of the text.

        Returns
        -------
        Results : `None`
            The text is added to the colorbar in place.
        """
        if color is None:
            color = "black"
        vmid = (cb.vmin + cb.vmax) / 2
        positions = {"vertical": (0.5, vmid), "horizontal": (vmid, 0.5)}
        cbtext = cb.ax.text(
            *positions[orientation],
            text,
            color=color,
            va="center",
            ha="center",
            fontsize=fontsize,
            fontweight=fontweight,
            rotation=orientation,
            alpha=alpha,
        )
        # Add a distinct outline around the text for better visibility in
        # various backgrounds.
        cbtext.set_path_effects(
            [mpl_path_effects.Stroke(linewidth=4, foreground="white", alpha=0.8), mpl_path_effects.Normal()]
        )

    @staticmethod
    def getLongestSuffixMatch(s, options):
        """Find the longest suffix in the provided list that matches the end of
        the given string.

        Parameters
        ----------
        s : `str`
            The target string for which we want to find a matching suffix.
        options : `list` [`str`]
            A list of potential suffix strings to match against the target
            string `s`.

        Returns
        -------
        Results : `str`
            The longest matching suffix from the `options` list. If no match is
            found, returns `None`.
        """
        return next((opt for opt in sorted(options, key=len, reverse=True) if s.endswith(opt)), None)
