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


from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import ChoiceField, Field
from lsst.pex.config.configurableActions import ConfigurableActionField
from lsst.utils.plotting import set_rubin_plotstyle
from matplotlib.figure import Figure

from ...actions.keyedData import CalcCompletenessHistogramAction
from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar
from .plotUtils import addPlotInfo

__all__ = ("CompletenessHist",)


class CompletenessHist(PlotAction):
    """Makes plots of completeness and purity."""

    action = ConfigurableActionField[CalcCompletenessHistogramAction](
        doc="Action to compute completeness/purity",
    )
    mag_ref_label = Field[str](doc="Label for the completeness x axis.", default="Reference magnitude")
    mag_target_label = Field[str](doc="Label for the purity x axis.", default="Measured magnitude")
    percentiles_style = ChoiceField[str](
        doc="Style and locations for completeness threshold percentile labels",
        allowed={
            "above_plot": "Labels in a semicolon-separated list above plot",
            "below_line": "Labels under the horizontal part of each line",
        },
        default="below_line",
    )
    publicationStyle = Field[bool](doc="Make a publication-style of plot", default=False)
    show_purity = Field[bool](doc="Whether to include a purity plot below completness", default=True)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.action.getOutputSchema()

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

    def makePlot(self, data, plotInfo, **kwargs):
        """Makes a plot showing the fraction of injected sources recovered by
        input magnitude.

        The behavior of this plot is controlled by `self.action`. This action
        must be added to a struct (usually self.process.calculateActions) by
        the tool that calls this plot.

        Parameters
        ----------
        data : `KeyedData`
            All the data
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:
            ``camera``
                The camera used to take the data (`lsst.afw.cameraGeom.Camera`)
            ``"cameraName"``
                The name of camera used to take the data (`str`).
            ``"filter"``
                The filter used for this data (`str`).
            ``"ccdKey"``
                The ccd/dectector key associated with this camera (`str`).
            ``"visit"``
                The visit of the data; only included if the data is from a
                single epoch dataset (`str`).
            ``"patch"``
                The patch that the data is from; only included if the data is
                from a coadd dataset (`str`).
            ``"tract"``
                The tract that the data comes from (`str`).
            ``"photoCalibDataset"``
                The dataset used for the calibration, e.g. "jointcal" or "fgcm"
                (`str`).
            ``"skyWcsDataset"``
                The sky Wcs dataset used (`str`).
            ``"rerun"``
                The rerun the data is stored in (`str`).

        Returns
        ------
        ``fig``
            The figure to be saved (`matplotlib.figure.Figure`).

        Notes
        -----
        The behaviour of this plot is largel

        Examples
        --------
        An example of the plot produced from this code from tract 3828 of the
        DC2 simulations is here:

        .. image:: /_static/analysis_tools/completenessPlotExample.png

        """

        # Make plot showing the fraction recovered in magnitude bins
        set_rubin_plotstyle()
        n_sub = 1 + self.show_purity
        fig, axes = plt.subplots(dpi=300, nrows=n_sub, figsize=(8, 4 * n_sub))
        if not self.show_purity:
            axes = [axes]
        color_counts = "purple"
        color_wrong = "firebrick"
        color_right = "teal"
        max_left = 1.05

        band = kwargs.get("band")
        action_hist = self.action.action
        names = {}
        for name in (
            "range_minimum",
            "range_maximum",
            "count",
            "count_ref",
            "count_target",
            "completeness",
            "completeness_bad_match",
            "completeness_good_match",
            "purity",
            "purity_bad_match",
            "purity_good_match",
        ):
            key = getattr(action_hist, f"name_{name}")
            if band is not None:
                key = key.format(band=band)
            names[name] = key

        ranges_min = data[names["range_minimum"]]
        ranges_max = data[names["range_maximum"]]
        x = (ranges_max + ranges_min) / 2.0
        interval = self.action.bins.mag_width / 1000.0
        x_err = interval / 2.0

        counts_all = data[names["count"]]

        plots = {
            "Completeness": {
                "count_type": "Reference",
                "counts": data[names["count_ref"]],
                "lines": (
                    (data[names["completeness"]], True, "k", "completeness"),
                    (data[names["completeness_bad_match"]], False, color_wrong, "wrong class"),
                    (data[names["completeness_good_match"]], False, color_right, "right class"),
                ),
                "xlabel": self.mag_ref_label,
            },
        }
        if self.show_purity:
            plots["Purity"] = {
                "count_type": "Object",
                "counts": data[names["count_target"]],
                "lines": (
                    (data[names["purity"]], True, "k", None),
                    (data[names["purity_bad_match"]], False, color_wrong, "wrong class"),
                    (data[names["purity_good_match"]], False, color_right, "right class"),
                ),
                "xlabel": self.mag_target_label,
            }

        # idx == 0 should be completeness; update this if that assumption
        # is changed
        for idx, (ylabel, plot_data) in enumerate(plots.items()):
            axes_idx = axes[idx]
            xlim = (ranges_min[0], ranges_max[-1])
            axes_idx.set(
                xlabel=plot_data["xlabel"],
                ylabel=ylabel,
                xlim=xlim,
                ylim=(0, max_left),
                xticks=np.arange(round(xlim[0]), round(xlim[1])),
                yticks=np.linspace(0, 1, 11),
            )
            axes_idx.grid(color="lightgrey", ls="-")
            ax_right = axes_idx.twinx()
            ax_right.set_ylabel(f"{plot_data['count_type']} counts/mag")
            ax_right.set_yscale("log")

            for y, do_err, color, label in plot_data["lines"]:
                axes_idx.errorbar(
                    x=x,
                    y=y,
                    xerr=x_err if do_err else None,
                    yerr=1.0 / np.sqrt(counts_all + 1) if do_err else None,
                    capsize=0,
                    color=color,
                    label=label,
                )
            y = plot_data["counts"] / interval
            # It should be unusual for np.max(y) to be zero; nonetheless...
            ax_right.step(
                [x[0] - interval] + list(x) + [x[-1] + interval],
                [0] + list(y) + [0],
                where="mid",
                color=color_counts,
                label="counts",
            )
            ax_right.set_ylim(0.999, 10 ** (max_left * np.log10(max(np.nanmax(y), 2))))
            ax_right.tick_params(axis="y", labelcolor=color_counts)
            lines_left, labels_left = axes_idx.get_legend_handles_labels()
            lines_right, labels_right = ax_right.get_legend_handles_labels()
            axes_idx.legend(lines_left + lines_right, labels_left + labels_right, loc="lower left", ncol=2)

            if idx == 0:
                percentiles = self.action.config_metrics.completeness_percentiles
                if percentiles:
                    above_plot = self.percentiles_style == "above_plot"
                    below_line = self.percentiles_style == "below_line"
                    kwargs_lines = dict(color="dimgrey", ls=":")
                    xlims = axes_idx.get_xlim()
                    if above_plot:
                        texts = []
                    elif below_line:
                        offset = 0.1 * (xlims[1] - xlims[0])
                    else:
                        raise RuntimeError(f"Unimplemented {self.percentiles_style=}")
                    for pct in percentiles:
                        name_pct = self.action.action.name_mag_completeness(
                            self.action.getPercentileName(pct),
                        )
                        if band is not None:
                            name_pct = name_pct.format(band=band)
                        mag_completeness = data.get(name_pct, None)
                        pct /= 100.0
                        if mag_completeness is not None and np.isfinite(mag_completeness):
                            axes_idx.plot([xlims[0], mag_completeness], [pct, pct], **kwargs_lines)
                            axes_idx.plot([mag_completeness, mag_completeness], [0, pct], **kwargs_lines)
                            text = f"{pct*100:.2g}%: {mag_completeness:.2f}"
                            if above_plot:
                                texts.append(text)
                            elif below_line:
                                axes_idx.text(
                                    mag_completeness - offset,
                                    pct,
                                    text,
                                    ha="right",
                                    va="top",
                                    fontsize=12,
                                )
                    if above_plot:
                        texts = f"Thresholds: {'; '.join(texts)}"
                        axes_idx.text(xlims[0], max_left, texts, ha="left", va="bottom")

        # Add useful information to the plot
        if not self.publicationStyle:
            addPlotInfo(fig, plotInfo)
        fig.tight_layout()
        fig.subplots_adjust(top=0.90)
        return fig
