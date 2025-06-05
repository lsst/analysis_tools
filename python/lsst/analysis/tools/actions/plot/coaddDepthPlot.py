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

__all__ = ("CoaddDepthPlot",)

from typing import TYPE_CHECKING, Any, Mapping
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from lsst.pex.config import ListField
from lsst.utils.plotting import publication_plots, set_rubin_plotstyle

from ...interfaces import PlotAction, Vector
from ..vector import BandSelector, PatchSelector
from .plotUtils import addPlotInfo

if TYPE_CHECKING:
    from ...interfaces import KeyedData, KeyedDataSchema

bands_dict = publication_plots.get_band_dicts()


class CoaddDepthPlot(PlotAction):
    """Make a plot of pixels per coadd depth."""

    threshold_list = ListField[int](
        default=[1, 3, 5, 12],
        doc="The n_image pixel value thresholds, in ascending order.",
    )

    def setDefaults(self):
        super().setDefaults()

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector]]] = []
        base.append(("patch", Vector))
        base.append(("band", Vector))
        base.append(("depth", Vector))
        base.append(("pixels", Vector))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Figure:
        self._validateInput(data)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData) -> None:
        needed = set(k[0] for k in self.getInputSchema())
        if not needed.issubset(data.keys()):
            raise ValueError(f"Input data does not contain all required keys: {self.getInputSchema()}")

    def makePlot(self, data: KeyedData, plotInfo: Mapping[str, str] | None = None, **kwargs: Any) -> Figure:
        """Make the plot.

        Parameters
        ----------
        `KeyedData`
            The catalog to plot the points from.

        plotInfo : `dict`
            A dictionary of the plot information.

        Returns
        -------
        fig : `~matplotlib.figure.Figure`
            The resulting figure.
        """
        set_rubin_plotstyle()
        fig = plt.figure(dpi=300, figsize=(20, 20))

        max_depth = max(data['depth'])
        max_pixels = max(data['pixels'])

        plt.subplots_adjust(hspace=0, wspace=0)

        patch_counter = 90  # The top left corner of a tract is patch 90
        m = 0  # subplot index
        while patch_counter >= 0:
            for n in range(10):  # column index
                ax = plt.subplot(10, 10, m + 1) # there are 10x10 patches per tract
                patchSelector = PatchSelector(vectorKey='patch', patches=[patch_counter])
                patch_mask = patchSelector(data)

                if patch_counter in data['patch']:
                    uniqueBands = set(data['band'][patch_mask])
                    for band in uniqueBands:
                        color = bands_dict['colors'][band]
                        markerstyle = bands_dict['symbols'][band]
                        bandSelector = BandSelector(vectorKey='band', bands=[band])
                        band_mask = bandSelector(data)

                        tot_mask = (patch_mask) & (band_mask)

                        ax.plot(data['depth'][tot_mask], data['pixels'][tot_mask],
                                color=color, linewidth=0, ls=None,
                                marker=markerstyle, ms=4, alpha=0.75,
                                label=f'{band}')
                        ax.grid(alpha=0.5)

                    # # add thresholds
                    # for threshold in self.threshold_list:
                    #     plt.axvline(x=threshold, color='m', alpha=0.75,
                    #                 ls='--', lw=1)

                # chart formatting
                # ax.set_xscale('symlog', base=2) # Need a solution for when max_depth is high, but not all patches/bands have quality coverage.
                ax.set_yscale('log')
                ax.set_xlim(0, max_depth + 5)
                ax.set_ylim(5, max_pixels)
                # ax.set_xticks(np.arange(0, max_depth, 20)) # Can we somehow generalize this?
                ax.tick_params(axis="both", which="minor")
                ax.tick_params(axis='both', which="both", top=False, right=False)

                # only label axes of the farmost left and bottom row of plots
                if (n != 0):
                    ax.set_yticklabels([])
                    ax.tick_params(axis='y', which='both', length=0)
                if (patch_counter not in range(10)):
                    ax.set_xticklabels([])
                    ax.tick_params(axis='x', which='both', length=0)

                ax.set_title(f"patch {patch_counter}", y=0.85)
                patch_counter += 1
                m += 1
            patch_counter -= 2*(n+1)
        fig.supxlabel('Number of input visits (n_image depth)', y=0.075)
        fig.supylabel('Count (pixels)', x=0.075)
        legend_elements = [
            Line2D([0], [0], color=bands_dict['colors']['u'], lw=0, marker=bands_dict['symbols']['u'], label='u'),
            Line2D([0], [0], color=bands_dict['colors']['g'], lw=0, marker=bands_dict['symbols']['g'], label='g'),
            Line2D([0], [0], color=bands_dict['colors']['r'], lw=0, marker=bands_dict['symbols']['r'], label='r'),
            Line2D([0], [0], color=bands_dict['colors']['i'], lw=0, marker=bands_dict['symbols']['i'], label='i'),
            Line2D([0], [0], color=bands_dict['colors']['z'], lw=0, marker=bands_dict['symbols']['z'], label='z'),
            Line2D([0], [0], color=bands_dict['colors']['y'], lw=0, marker=bands_dict['symbols']['y'], label='y'),
        ]
        plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 10))

        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
