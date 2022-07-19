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

__all__ = ("HistPanel", "HistPlot")

from collections import defaultdict
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import Config, ConfigDictField, DictField, Field
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from scipy.stats import median_absolute_deviation as sigmaMad

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from .plotUtils import addPlotInfo


class HistPanel(Config):
    label = Field[str](
        doc="Panel x-axis label.",
        default="label",
    )
    hists = DictField[str, str](
        doc="A dict specifying the histograms to be plotted in this panel. Keys are used to identify "
        "histogram IDs. Values are used to add to the legend label displayed in the upper corner of the "
        "panel.",
        optional=False,
    )
    yscale = Field[str](
        doc="Y axis scaling.",
        default="linear",
    )
    bins = Field[int](
        doc="Number of x axis bins.",
        default=50,
    )
    pLower = Field[float](
        doc="Percentile used to determine the lower range of the histogram bins. If more than one histogram "
        "is plotted in the panel, the percentile limit is the minimum value across all input data.",
        default=2.0,
    )
    pUpper = Field[float](
        doc="Percentile used to determine the upper range of the histogram bins. If more than one histogram "
        "is plotted, the percentile limit is the maximum value across all input data.",
        default=98.0,
    )


class HistPlot(PlotAction):
    panels = ConfigDictField(
        doc="A configurable dict describing the panels to be plotted, and the histograms for each panel.",
        keytype=str,
        itemtype=HistPanel,
        default={},
    )
    cmap = Field[str](
        doc="Color map used for histogram lines. All types available via `plt.cm` may be used. "
        "A number of custom color maps are also defined: `newtab10`, `bright`, `vibrant`.",
        default="newtab10",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        for panel in self.panels:  # type: ignore
            for histData in self.panels[panel].hists.items():  # type: ignore
                yield histData, Vector

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        return self.makePlot(data, **kwargs)
        # table is a dict that needs: x, y, run, skymap, filter, tract,

    def makePlot(
        self, data: KeyedData, plotInfo: Mapping[str, str] = None, **kwargs  # type: ignore
    ) -> Figure:
        """Make an N-panel plot with a user-configurable number of histograms
        displayed in each panel.

        Parameters
        ----------
        data : `pandas.core.frame.DataFrame`
            The catalog to plot the points from.
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:
                `"run"`
                    Output run for the plots (`str`).
                `"tractTableType"`
                    Table from which results are taken (`str`).
                `"plotName"`
                    Output plot name (`str`)
                `"SN"`
                    The global signal-to-noise data threshold (`float`)
                `"skymap"`
                    The type of skymap used for the data (`str`).
                `"tract"`
                    The tract that the data comes from (`int`).
                `"bands"`
                    The bands used for this data (`str` or `list`).
                `"visit"`
                    The visit that the data comes from (`int`)

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.

        """

        # set up figure
        fig = plt.figure(dpi=300)
        hist_fig, side_fig = fig.subfigures(1, 2, wspace=0, width_ratios=[3, 1])
        axs = self._makeAxes(hist_fig)

        # loop over each panel; plot histograms
        cols = self._assignColors()
        all_handles, all_nums, all_meds, all_mads = [], [], [], []
        for panel, ax in zip(self.panels, axs):
            nums, meds, mads = self._makePanel(data, panel, ax, cols[panel], **kwargs)
            handles, labels = ax.get_legend_handles_labels()  # code for plotting
            all_handles += handles
            all_nums += nums
            all_meds += meds
            all_mads += mads

        # add side panel; add statistics
        self._addStatisticsPanel(side_fig, all_handles, all_nums, all_meds, all_mads)

        # add general plot info
        hist_fig = addPlotInfo(hist_fig, plotInfo)

        # finish up
        hist_fig.text(0.01, 0.42, "Frequency", rotation=90, transform=hist_fig.transFigure)
        plt.draw()
        return fig

    def _makeAxes(self, fig):
        """Determine axes layout for main histogram figure."""
        num_panels = len(self.panels)
        if num_panels <= 1:
            ncols = 1
        else:
            ncols = 2
        nrows = int(np.ceil(num_panels / ncols))

        gs = GridSpec(nrows, ncols, left=0.13, right=0.99, bottom=0.1, top=0.88, wspace=0.25, hspace=0.45)

        axs = []
        counter = 0
        for row in range(nrows):
            for col in range(ncols):
                counter += 1
                if counter < num_panels:
                    axs.append(fig.add_subplot(gs[row : row + 1, col : col + 1]))
                else:
                    axs.append(fig.add_subplot(gs[row : row + 1, col : np.min([col + 2, ncols + 1])]))
                    break

        return axs

    def _assignColors(self):
        """Assign colors to histograms using a given color map."""
        custom_cmaps = dict(
            # https://www.tableau.com/about/blog/2016/7/colors-upgrade-tableau-10-56782
            newtab10=[
                "#4e79a7",
                "#f28e2b",
                "#e15759",
                "#76b7b2",
                "#59a14f",
                "#edc948",
                "#b07aa1",
                "#ff9da7",
                "#9c755f",
                "#bab0ac",
            ],
            # https://personal.sron.nl/~pault/#fig:scheme_bright
            bright=[
                "#4477AA",
                "#EE6677",
                "#228833",
                "#CCBB44",
                "#66CCEE",
                "#AA3377",
                "#BBBBBB",
            ],
            # https://personal.sron.nl/~pault/#fig:scheme_vibrant
            vibrant=[
                "#EE7733",
                "#0077BB",
                "#33BBEE",
                "#EE3377",
                "#CC3311",
                "#009988",
                "#BBBBBB",
            ],
        )
        if self.cmap in custom_cmaps.keys():
            all_cols = custom_cmaps[self.cmap]
        else:
            try:
                all_cols = getattr(plt.cm, self.cmap).copy().colors
            except AttributeError:
                raise ValueError(f"Unrecognized color map: {self.cmap}")

        counter = 0
        cols = defaultdict(list)
        for panel in self.panels:
            for hist in self.panels[panel].hists:
                cols[panel].append(all_cols[counter % len(all_cols)])
                counter += 1
        return cols

    def _makePanel(self, data, panel, ax, col, **kwargs):
        """Plot a single panel containing histograms."""
        panel_range = self._getPanelRange(data, panel)
        nums, meds, mads = [], [], []
        for i, hist in enumerate(self.panels[panel].hists):
            hist_data = data[hist][np.isfinite(data[hist])]
            ax.hist(
                hist_data,
                range=panel_range,
                bins=self.panels[panel].bins,
                histtype="step",
                lw=2,
                color=col[i],
                label=self.panels[panel].hists[hist],
            )
            num, med, mad = self._calcStats(hist_data)
            nums.append(num)
            meds.append(med)
            mads.append(mad)
            ax.axvline(med, ls="--", lw=1, c=col[i])
        ax.legend(fontsize=6, loc="upper left")
        ax.set_xlim(panel_range)
        ax.set_xlabel(self.panels[panel].label)
        ax.set_yscale(self.panels[panel].yscale)
        ax.tick_params(labelsize=7)
        # add a buffer to the top of the plot to allow headspace for labels
        ylims = list(ax.get_ylim())
        if ax.get_yscale() == "log":
            ylims[1] = 10 ** (np.log10(ylims[1]) * 1.1)
        else:
            ylims[1] *= 1.1
        ax.set_ylim(ylims[0], ylims[1])
        return nums, meds, mads

    def _getPanelRange(self, data, panel):
        """Determine panel x-axis range based on data percentile limits."""
        panel_range = [np.nan, np.nan]
        for hist in self.panels[panel].hists:
            hist_range = np.nanpercentile(data[hist], [self.panels[panel].pLower, self.panels[panel].pUpper])
            panel_range[0] = np.nanmin([panel_range[0], hist_range[0]])
            panel_range[1] = np.nanmax([panel_range[1], hist_range[1]])
        return panel_range

    def _calcStats(self, data):
        """Calculate the number of data points, median, and median absolute
        deviation of input data."""
        num = len(data)
        med = np.nanmedian(data)
        mad = sigmaMad(data)
        return num, med, mad

    def _addStatisticsPanel(self, fig, handles, nums, meds, mads):
        """Add an adjoining panel containing histogram summary statistics."""
        ax = fig.add_subplot(1, 1, 1)
        ax.axis("off")
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.09, top=0.9)

        # empty handle, used to populate the bespoke legend layout
        empty = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0)

        # set up new legend handles and labels
        legend_handles = [empty] + handles + ([empty] * 3 * len(handles)) + ([empty] * 3)
        legend_labels = (
            ([""] * (len(handles) + 1))
            + ["Num"]
            + nums
            + ["Med"]
            + [f"{x:0.1f}" for x in meds]
            + ["${{\\sigma}}_{{MAD}}$"]
            + [f"{x:0.1f}" for x in mads]
        )

        # add the legend
        ax.legend(
            legend_handles,
            legend_labels,
            loc="lower left",
            ncol=4,
            handletextpad=-0.25,
            fontsize=6,
            borderpad=0,
            frameon=False,
            columnspacing=-0.25,
        )
