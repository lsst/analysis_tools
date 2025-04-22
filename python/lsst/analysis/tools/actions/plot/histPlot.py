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

__all__ = ("HistPanel", "HistPlot", "HistStatsPanel")

import logging
from collections import defaultdict
from typing import Mapping

import numpy as np
from lsst.pex.config import (
    ChoiceField,
    Config,
    ConfigDictField,
    ConfigField,
    DictField,
    Field,
    FieldValidationError,
    ListField,
)
from lsst.utils.plotting import make_figure
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from ...math import nanMax, nanMedian, nanMin, sigmaMad
from .plotUtils import addPlotInfo

log = logging.getLogger(__name__)


class HistStatsPanel(Config):
    """A Config class that holds parameters to configure a the stats panel
    shown for histPlot.

    The fields in this class correspond to the parameters that can be used to
    customize the HistPlot stats panel.

    - The ListField parameter a dict to specify names of 3 stat columns accepts
      latex formating

    - The other parameters (stat1, stat2, stat3) are lists of strings that
      specify vector keys correspoinding to scalar values computed in the
      prep/process/produce steps of an analysis tools plot/metric configurable
      action. There should be one key for each group in the HistPanel.

    A separate config class is used instead of constructing
    `~lsst.pex.config.DictField`'s in HistPanel for each parameter for clarity
    and consistency.

    Notes
    -----
    This is intended to be used as a configuration of the HistPlot/HistPanel
    class.

    If no HistStatsPanel is specified then the default behavor persists where
    the stats panel shows N / median / sigma_mad for each group in the panel.
    """

    statsLabels = ListField[str](
        doc="list specifying the labels for stats",
        length=3,
        default=("N$_{{data}}$", "Med", "${{\\sigma}}_{{MAD}}$"),
    )
    stat1 = ListField[str](
        doc="A list specifying the vector keys of the first scalar statistic to be shown in this panel."
        "there should be one entry for each hist in the panel",
        default=None,
        optional=True,
    )
    stat2 = ListField[str](
        doc="A list specifying the vector keys of the second scalar statistic to be shown in this panel."
        "there should be one entry for each hist in the panel",
        default=None,
        optional=True,
    )
    stat3 = ListField[str](
        doc="A list specifying the vector keys of the third scalar statistic to be shown in this panel."
        "there should be one entry for each hist in the panel",
        default=None,
        optional=True,
    )

    def validate(self):
        super().validate()
        if not all([self.stat1, self.stat2, self.stat3]) and any([self.stat1, self.stat2, self.stat3]):
            raise ValueError(f"{self._name}: If one stat is configured, all 3 stats must be configured")


class HistPanel(Config):
    """A Config class that holds parameters to configure a single panel of a
    histogram plot. This class is intended to be used within the ``HistPlot``
    class.
    """

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
        doc="Number of x axis bins within plot x-range.",
        default=50,
    )
    rangeType = ChoiceField[str](
        doc="Set the type of range to use for the x-axis. Range bounds will be set according to "
        "the values of lowerRange and upperRange.",
        allowed={
            "percentile": "Upper and lower percentile ranges of the data.",
            "sigmaMad": "Range is (sigmaMad - lowerRange*sigmaMad, sigmaMad + upperRange*sigmaMad).",
            "fixed": "Range is fixed to (lowerRange, upperRange).",
        },
        default="percentile",
    )
    lowerRange = Field[float](
        doc="Lower range specifier for the histogram bins. See rangeType for interpretation "
        "based on the type of range requested. If more than one histogram is plotted in a given "
        "panel and rangeType is not set to fixed, the limit is the minimum value across all input "
        "data.",
        default=0.0,
    )
    upperRange = Field[float](
        doc="Upper range specifier for the histogram bins. See rangeType for interpretation "
        "based on the type of range requested. If more than one histogram is plotted in a given "
        "panel and rangeType is not set to fixed, the limit is the maximum value across all input "
        "data.",
        default=100.0,
    )
    referenceValue = Field[float](
        doc="Value at which to add a black solid vertical line. Ignored if set to `None`.",
        default=None,
        optional=True,
    )
    refRelativeToMedian = Field[bool](
        doc="Is the referenceValue meant to be an offset from the median?",
        default=False,
        optional=True,
    )
    histDensity = Field[bool](
        doc="Whether to plot the histogram as a normalized probability distribution. Must also "
        "provide a value for referenceValue",
        default=False,
    )
    statsPanel = ConfigField[HistStatsPanel](
        doc="configuration for stats to be shown on plot, if None then "
        "default stats: N, median, sigma mad are shown",
        default=None,
    )

    def validate(self):
        super().validate()
        if self.rangeType == "percentile" and self.lowerRange < 0.0 or self.upperRange > 100.0:
            msg = (
                "For rangeType %s, ranges must obey: lowerRange >= 0 and upperRange <= 100." % self.rangeType
            )
            raise FieldValidationError(self.__class__.rangeType, self, msg)
        if self.rangeType == "sigmaMad" and self.lowerRange < 0.0:
            msg = (
                "For rangeType %s, lower range must obey: lowerRange >= 0 (the lower range is "
                "set as median - lowerRange*sigmaMad." % self.rangeType
            )
            raise FieldValidationError(self.__class__.rangeType, self, msg)
        if self.rangeType == "fixed" and (self.upperRange - self.lowerRange) == 0.0:
            msg = (
                "For rangeType %s, lower and upper ranges must differ (i.e. must obey: "
                "upperRange - lowerRange != 0)." % self.rangeType
            )
            raise FieldValidationError(self.__class__.rangeType, self, msg)
        if self.histDensity and self.referenceValue is None:
            msg = "Must provide referenceValue if histDensity is True."
            raise FieldValidationError(self.__class__.referenceValue, self, msg)


class HistPlot(PlotAction):
    """Make an N-panel plot with a configurable number of histograms displayed
    in each panel. Reference lines showing values of interest may also be added
    to each histogram. Panels are configured using the ``HistPanel`` class.
    """

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

        Examples
        --------
        An example histogram plot may be seen below:

        .. image:: /_static/analysis_tools/histPlotExample.png

        For further details on how to generate a plot, please refer to the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """

        # set up figure
        fig = make_figure(dpi=300)
        hist_fig, side_fig = fig.subfigures(1, 2, wspace=0, width_ratios=[2.9, 1.1])
        axs, ncols, nrows = self._makeAxes(hist_fig)

        # loop over each panel; plot histograms
        colors = self._assignColors()
        nth_panel = len(self.panels)
        nth_col = ncols
        nth_row = nrows - 1
        label_font_size = max(6, 10 - nrows)
        for panel, ax in zip(self.panels, axs):
            nth_panel -= 1
            nth_col = ncols - 1 if nth_col == 0 else nth_col - 1
            if nth_panel == 0 and nrows * ncols - len(self.panels) > 0:
                nth_col -= 1
            # Set font size for legend based on number of panels being plotted.
            legend_font_size = max(4, int(7 - len(self.panels[panel].hists) / 2 - nrows // 2))  # type: ignore
            nums, meds, mads, stats_dict = self._makePanel(
                data,
                panel,
                ax,
                colors[panel],
                label_font_size=label_font_size,
                legend_font_size=legend_font_size,
                ncols=ncols,
            )

            all_handles, all_nums, all_meds, all_mads = [], [], [], []
            handles, labels = ax.get_legend_handles_labels()  # code for plotting
            all_handles += handles
            all_nums += nums
            all_meds += meds
            all_mads += mads
            title_str = self.panels[panel].label  # type: ignore
            # add side panel; add statistics
            self._addStatisticsPanel(
                side_fig,
                all_handles,
                all_nums,
                all_meds,
                all_mads,
                stats_dict,
                legend_font_size=legend_font_size,
                yAnchor0=ax.get_position().y0,
                nth_row=nth_row,
                nth_col=nth_col,
                title_str=title_str,
            )
            nth_row = nth_row - 1 if nth_col == 0 else nth_row

        # add general plot info
        if plotInfo is not None:
            hist_fig = addPlotInfo(hist_fig, plotInfo)

        # finish up
        fig.canvas.draw()
        return fig

    def _makeAxes(self, fig):
        """Determine axes layout for main histogram figure."""
        num_panels = len(self.panels)
        if num_panels <= 1:
            ncols = 1
        else:
            ncols = 2
        nrows = int(np.ceil(num_panels / ncols))

        gs = GridSpec(nrows, ncols, left=0.12, right=0.88, bottom=0.1, top=0.88, wspace=0.41, hspace=0.45)

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

        return axs, ncols, nrows

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
            all_colors = custom_cmaps[self.cmap]
        else:
            try:
                all_colors = getattr(cm, self.cmap).copy().colors
            except AttributeError:
                raise ValueError(f"Unrecognized color map: {self.cmap}")

        counter = 0
        colors = defaultdict(list)
        for panel in self.panels:
            for hist in self.panels[panel].hists:
                colors[panel].append(all_colors[counter % len(all_colors)])
                counter += 1
        return colors

    def _makePanel(self, data, panel, ax, colors, label_font_size=9, legend_font_size=7, ncols=1):
        """Plot a single panel containing histograms."""
        nums, meds, mads = [], [], []
        for i, hist in enumerate(self.panels[panel].hists):
            hist_data = data[hist][np.isfinite(data[hist])]
            num, med, mad = self._calcStats(hist_data)
            nums.append(num)
            meds.append(med)
            mads.append(mad)
        panel_range = self._getPanelRange(data, panel, mads=mads, meds=meds)
        if all(np.isfinite(panel_range)):
            nHist = 0
            for i, hist in enumerate(self.panels[panel].hists):
                hist_data = data[hist][np.isfinite(data[hist])]
                if len(hist_data) > 0:
                    ax.hist(
                        hist_data,
                        range=panel_range,
                        bins=self.panels[panel].bins,
                        histtype="step",
                        density=self.panels[panel].histDensity,
                        lw=2,
                        color=colors[i],
                        label=self.panels[panel].hists[hist],
                    )
                    ax.axvline(meds[i], ls=(0, (5, 3)), lw=1, c=colors[i])
                    nHist += 1

            if nHist > 0:
                ax.legend(fontsize=legend_font_size, loc="upper left", frameon=False)
            ax.set_xlim(panel_range)
            # The following accommodates spacing for ranges with large numbers
            # but small-ish dynamic range (example use case: RA 300-301).
            if ncols > 1 and max(np.abs(panel_range)) >= 100 and (panel_range[1] - panel_range[0]) < 5:
                ax.xaxis.set_major_formatter("{x:.2f}")
                ax.tick_params(axis="x", labelrotation=25, pad=-1)
            ax.set_xlabel(self.panels[panel].label, fontsize=label_font_size)
            y_label = "Normalized (PDF)" if self.panels[panel].histDensity else "Frequency"
            ax.set_ylabel(y_label, fontsize=label_font_size)
            ax.set_yscale(self.panels[panel].yscale)
            ax.tick_params(labelsize=max(5, label_font_size - 2))
            # add a buffer to the top of the plot to allow headspace for labels
            ylims = list(ax.get_ylim())
            if ax.get_yscale() == "log":
                ylims[1] = 10 ** (np.log10(ylims[1]) * 1.1)
            else:
                ylims[1] *= 1.1
            ax.set_ylim(ylims[0], ylims[1])

            # Draw a vertical line at a reference value, if given.
            # If histDensity is True, also plot a reference PDF with
            # mean = referenceValue and sigma = 1 for reference.
            if self.panels[panel].referenceValue is not None:
                ax = self._addReferenceLines(ax, panel, panel_range, meds, legend_font_size=legend_font_size)

            # Check if we should use the default stats panel or if a custom one
            # has been created.
            statList = [
                self.panels[panel].statsPanel.stat1,
                self.panels[panel].statsPanel.stat2,
                self.panels[panel].statsPanel.stat3,
            ]
            if not any(statList):
                stats_dict = {
                    "statLabels": ["N$_{{data}}$", "Med", "${{\\sigma}}_{{MAD}}$"],
                    "stat1": nums,
                    "stat2": meds,
                    "stat3": mads,
                }
            elif all(statList):
                stat1 = [data[stat] for stat in self.panels[panel].statsPanel.stat1]
                stat2 = [data[stat] for stat in self.panels[panel].statsPanel.stat2]
                stat3 = [data[stat] for stat in self.panels[panel].statsPanel.stat3]
                stats_dict = {
                    "statLabels": self.panels[panel].statsPanel.statsLabels,
                    "stat1": stat1,
                    "stat2": stat2,
                    "stat3": stat3,
                }
            else:
                raise RuntimeError("Invalid configuration of HistStatPanel")
        else:
            stats_dict = {key: [] for key in ("stat1", "stat2", "stat3")}
            stats_dict["statLabels"] = [""] * 3
        return nums, meds, mads, stats_dict

    def _getPanelRange(self, data, panel, mads=None, meds=None):
        """Determine panel x-axis range based config settings."""
        panel_range = [np.nan, np.nan]
        rangeType = self.panels[panel].rangeType
        lowerRange = self.panels[panel].lowerRange
        upperRange = self.panels[panel].upperRange
        if rangeType == "percentile":
            panel_range = self._getPercentilePanelRange(data, panel)
        elif rangeType == "sigmaMad":
            # Set the panel range to extend lowerRange[upperRange] times the
            # maximum sigmaMad for the datasets in the panel to the left[right]
            # from the minimum[maximum] median value of all datasets in the
            # panel.
            maxMad = nanMax(mads)
            maxMed = nanMax(meds)
            minMed = nanMin(meds)
            panel_range = [minMed - lowerRange * maxMad, maxMed + upperRange * maxMad]
            if panel_range[1] - panel_range[0] == 0:
                log.info(
                    "NOTE: panel_range for {} based on med/sigMad was 0. Computing using "
                    "percentile range instead.".format(panel)
                )
                panel_range = self._getPercentilePanelRange(data, panel)
        elif rangeType == "fixed":
            panel_range = [lowerRange, upperRange]
        else:
            raise RuntimeError(f"Invalid rangeType: {rangeType}")
        return panel_range

    def _getPercentilePanelRange(self, data, panel):
        """Determine panel x-axis range based on data percentile limits."""
        panel_range = [np.nan, np.nan]
        for hist in self.panels[panel].hists:
            data_hist = data[hist]
            # TODO: Consider raising instead
            if len(data_hist) > 0:
                hist_range = np.nanpercentile(
                    data[hist], [self.panels[panel].lowerRange, self.panels[panel].upperRange]
                )
                panel_range[0] = nanMin([panel_range[0], hist_range[0]])
                panel_range[1] = nanMax([panel_range[1], hist_range[1]])
        return panel_range

    def _calcStats(self, data):
        """Calculate the number of data points, median, and median absolute
        deviation of input data."""
        num = len(data)
        med = nanMedian(data)
        mad = sigmaMad(data)
        return num, med, mad

    def _addReferenceLines(self, ax, panel, panel_range, meds, legend_font_size=7):
        """Draw the vertical reference line and density curve (if requested)
        on the panel.
        """
        ax2 = ax.twinx()
        ax2.axis("off")
        ax2.set_xlim(ax.get_xlim())
        ax2.set_ylim(ax.get_ylim())

        if self.panels[panel].histDensity:
            reference_label = None
        else:
            if self.panels[panel].refRelativeToMedian:
                reference_value = self.panels[panel].referenceValue + meds[0]
                reference_label = "${{\\mu_{{ref}}}}$: {:10.3F}".format(reference_value)
            else:
                reference_value = self.panels[panel].referenceValue
                reference_label = "${{\\mu_{{ref}}}}$: {:10.3F}".format(reference_value)
            ax2.axvline(reference_value, ls="-", lw=1, c="black", zorder=0, label=reference_label)
        if self.panels[panel].histDensity:
            ref_x = np.arange(panel_range[0], panel_range[1], (panel_range[1] - panel_range[0]) / 100.0)
            ref_mean = self.panels[panel].referenceValue
            ref_std = 1.0
            ref_y = (
                1.0
                / (ref_std * np.sqrt(2.0 * np.pi))
                * np.exp(-((ref_x - ref_mean) ** 2) / (2.0 * ref_std**2))
            )
            ax2.fill_between(ref_x, ref_y, alpha=0.1, color="black", label="P$_{{norm}}(0,1)$", zorder=-1)
            # Make sure the y-axis extends beyond the data plotted and that
            # the y-ranges of both axes are in sync.
            y_max = max(max(ref_y), ax2.get_ylim()[1])
            if ax2.get_ylim()[1] < 1.05 * y_max:
                ax.set_ylim(ax.get_ylim()[0], 1.05 * y_max)
                ax2.set_ylim(ax.get_ylim())
        ax2.legend(fontsize=legend_font_size, handlelength=1.5, loc="upper right", frameon=False)

        return ax

    def _addStatisticsPanel(
        self,
        fig,
        handles,
        nums,
        meds,
        mads,
        stats_dict,
        legend_font_size=8,
        yAnchor0=0.0,
        nth_row=0,
        nth_col=0,
        title_str=None,
    ):
        """Add an adjoining panel containing histogram summary statistics."""
        ax = fig.add_subplot(1, 1, 1)
        ax.axis("off")
        fig.subplots_adjust(left=0.05, right=0.95, bottom=0.0, top=1.0)
        # empty handle, used to populate the bespoke legend layout
        empty = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0)

        # set up new legend handles and labels
        legend_handles = [empty] + handles + ([empty] * 3 * len(handles)) + ([empty] * 3)

        legend_labels = (
            ([""] * (len(handles) + 1))
            + [stats_dict["statLabels"][0]]
            + [f"{x:.3g}" if abs(x) > 0.01 else f"{x:.2e}" for x in stats_dict["stat1"]]
            + [stats_dict["statLabels"][1]]
            + [f"{x:.3g}" if abs(x) > 0.01 else f"{x:.2e}" for x in stats_dict["stat2"]]
            + [stats_dict["statLabels"][2]]
            + [f"{x:.3g}" if abs(x) > 0.01 else f"{x:.2e}" for x in stats_dict["stat3"]]
        )
        # Replace "e+0" with "e" and "e-0" with "e-" to save space.
        legend_labels = [label.replace("e+0", "e") for label in legend_labels]
        legend_labels = [label.replace("e-0", "e-") for label in legend_labels]

        # Set the y anchor for the legend such that it roughly lines up with
        # the panels.
        yAnchor = max(0, yAnchor0 - 0.01) + nth_col * (0.008 + len(nums) * 0.005) * legend_font_size

        nth_legend = ax.legend(
            legend_handles,
            legend_labels,
            loc="lower left",
            bbox_to_anchor=(-0.25, yAnchor),
            ncol=4,
            handletextpad=-0.25,
            fontsize=legend_font_size,
            borderpad=0,
            frameon=False,
            columnspacing=-0.25,
            title=title_str,
            title_fontproperties={"weight": "bold", "size": legend_font_size},
        )
        if nth_row + nth_col > 0:
            ax.add_artist(nth_legend)
