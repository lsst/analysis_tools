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

__all__ = ("BarPanel", "BarPlot")

import operator as op
from collections import defaultdict
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import Config, ConfigDictField, DictField, Field
from lsst.utils.plotting import set_rubin_plotstyle
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from .plotUtils import addPlotInfo


class BarPanel(Config):
    """A configurable class describing a panel in a bar plot."""

    label = Field[str](
        doc="Panel x-axis label.",
        default="label",
    )
    bars = DictField[str, str](
        doc="A dict specifying the bar graphs to be plotted in this panel. Keys are used to identify "
        "bar graph IDs. Values are used to add to the legend label displayed in the upper corner of the "
        "panel.",
        optional=False,
    )
    yscale = Field[str](
        doc="Y axis scaling.",
        default="linear",
    )


class BarPlot(PlotAction):
    """A plotting tool which can take multiple keyed data inputs
    and can create one or more bar graphs.
    """

    panels = ConfigDictField(
        doc="A configurable dict describing the panels to be plotted, and the bar graphs for each panel.",
        keytype=str,
        itemtype=BarPanel,
        default={},
    )
    cmap = Field[str](
        doc="Color map used for bar lines. All types available via `plt.cm` may be used. "
        "A number of custom color maps are also defined: `newtab10`, `bright`, `vibrant`.",
        default="newtab10",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        for panel in self.panels:  # type: ignore
            for barData in self.panels[panel].bars.items():  # type: ignore
                yield barData, Vector

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        return self.makePlot(data, **kwargs)

    def makePlot(
        self, data: KeyedData, plotInfo: Mapping[str, str] = None, **kwargs  # type: ignore
    ) -> Figure:
        """Make an N-panel plot with a user-configurable number of bar graphs
        displayed in each panel.

        Parameters
        ----------
        data : `KeyedData`
            The catalog to plot the points from.
        plotInfo : `dict`
            An optional dictionary of information about the data being
            plotted with keys:

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
        fig = plt.figure(dpi=400)
        bar_fig, side_fig = fig.subfigures(1, 2, wspace=0, width_ratios=[3, 1])
        axs = self._makeAxes(bar_fig)

        # Set the rubin style
        set_rubin_plotsyle()

        # loop over each panel; plot bar graphs
        cols = self._assignColors()
        all_handles, all_nums, all_vector_labels, all_x_values = [], [], [], []
        for panel, ax in zip(self.panels, axs):
            nums, sorted_label, sorted_x_values = self._makePanel(data, panel, ax, cols[panel], **kwargs)
            handles, labels = ax.get_legend_handles_labels()  # code for plotting
            all_handles += handles
            all_nums += nums
            all_vector_labels += sorted_label
            all_x_values += sorted_x_values

        # add side panel; add statistics
        self._addStatisticsPanel(side_fig, all_handles, all_nums, all_vector_labels, all_x_values)

        # add general plot info
        if plotInfo is not None:
            bar_fig = addPlotInfo(bar_fig, plotInfo)

        # finish up
        bar_fig.text(0.01, 0.42, "Frequency", rotation=90, transform=bar_fig.transFigure)
        plt.draw()
        return fig

    def _makeAxes(self, fig):
        """Determine axes layout for main bar graph figure."""
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
        """Assign colors to bar graphs using a given color map."""
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
            for bar in self.panels[panel].bars:
                cols[panel].append(all_cols[counter % len(all_cols)])
                counter += 1
        return cols

    def _makePanel(self, data, panel, ax, col, **kwargs):
        """Plot a single panel containing bar graphs."""
        nums = []
        x_values, assigned_labels, assigned_colors = self._assignBinElements(data, panel, col)
        sorted_x_values, sorted_labels, sorted_colors = self._sortBarBins(
            x_values, assigned_labels, assigned_colors
        )
        width, columns = self._getBarWidths(sorted_x_values)

        for i, bin in enumerate(sorted_x_values):
            bar_data = op.countOf(data[sorted_labels[i]][np.isfinite(data[sorted_labels[i]])], bin)

            if width[i] == 1:
                bin_center = bin
            else:
                bin_center = bin - 0.35 + width[i] * columns[i]

            ax.bar(bin_center, bar_data, width[i], lw=2, label=sorted_labels[i], color=sorted_colors[i])
            nums.append(bar_data)

        # Get plot range
        x_range = [x for x in range(int(min(sorted_x_values)), int(max(sorted_x_values)) + 1)]
        ax.set_xticks(x_range)
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
        return nums, sorted_labels, sorted_x_values

    def _assignBinElements(self, data, panel, col):
        labels = []
        assigned_labels = []
        x_values = []
        assigned_colors = []
        n_labels = 0

        for bar in self.panels[panel].bars:
            labels.append(bar)
            n_labels += 1

        # If a label has multiple unique elements in it, repeats the label
        i = 0
        for single_label in labels:
            unique_elements = np.unique(data[single_label])

            for bin in unique_elements:
                x_values.append(int(bin))

            for count in range(len(unique_elements)):
                assigned_labels.append(single_label)
                assigned_colors.append(col[i])  # Assign color from color cmap

            i += 1

        return x_values, assigned_labels, assigned_colors

    def _sortBarBins(self, x_values, assigned_labels, assigned_colors):
        """Sorts the existing x_values, assigned_labels,
        and assigned_colors/x_value from lowest to
        highest and then uses the sorted indices to sort
        all x, labels, and colors in that order.
        """

        sorted_indices = np.argsort(x_values)

        sorted_labels = []
        sorted_x_values = []
        sorted_colors = []

        for position in sorted_indices:
            sorted_x_values.append(x_values[position])
            sorted_labels.append(assigned_labels[position])
            sorted_colors.append(assigned_colors[position])

        return sorted_x_values, sorted_labels, sorted_colors

    def _getBarWidths(self, x_values):
        """Determine the width of the panels in each
        bin and which column is assigned."""
        width = []
        columns = []
        current_column = 0
        current_i = 0

        for i in x_values:
            # Number of repeating values
            n_repeating = x_values.count(i)
            width.append(1.0 / n_repeating)
            if n_repeating > 1 and current_column != 0 and current_i == i:
                columns.append(current_column)
                current_column += 1

            else:
                current_column = 0
                columns.append(current_column)
                current_i = i
                current_column += 1

        return width, columns

    def _addStatisticsPanel(self, fig, handles, nums, sorted_labels, sorted_x_value):
        """Add an adjoining panel containing bar graph summary statistics."""
        ax = fig.add_subplot(1, 1, 1)
        ax.axis("off")
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.09, top=0.9)

        # empty handle, used to populate the bespoke legend layout
        empty = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0)

        # set up new legend handles and labels

        legend_handles = [empty] + handles + ([empty] * 3 * len(handles)) + ([empty] * 3)
        legend_labels = (
            ([""] * (len(handles) + 1))
            + ["Bin"]
            + sorted_x_value
            + ["Count"]
            + nums
            + ["Sources"]
            + sorted_labels
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
