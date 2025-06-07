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

from lsst.utils.plotting import publication_plots

from ...interfaces import PlotAction, Vector
from ..vector import BandSelector, PatchSelector, UniqueAction
from .plotUtils import addPlotInfo

if TYPE_CHECKING:
    from ...interfaces import KeyedData, KeyedDataSchema

bands_dict = publication_plots.get_band_dicts()


class CoaddDepthPlot(PlotAction):
    """Make a plot (with quantiles) of pixels per coadd depth."""

    # quantile_list = ListField[str](
    #     doc="The percentiles at which to compute n_image values, in ascending order.",
    #     default=[5, 10, 25, 50, 75, 90, 95],
    #     dtype=int,
    # )

    def setDefaults(self):
        super().setDefaults()

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector]]] = []
        base.append(("patch", Vector))
        base.append(("band", Vector))
        base.append(("depth", Vector))
        base.append(("pixels", Vector))

        # base.append(("q_band", Vector))
        # base.append(("q_depth", Vector))
        # base.append(("quantiles", Vector))
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
        fig = plt.figure(dpi=300, figsize=(20, 20))

        max_depth = max(data['depth'])
        max_pixels = max(data['pixels'])

        plt.subplots_adjust(hspace=0, wspace=0)

        patch_counter = 90  # The top left corner of a tract is patch 90
        m = 0  # subplot index
        while patch_counter >= 0:
            for n in range(10):  # column index
                ax = plt.subplot(10, 10, m + 1)
                patchSelector = PatchSelector(vectorKey='patch', patches=[patch_counter])
                patch_mask = patchSelector(data)

                # uniqueBands = UniqueAction(vectorKey='band')  # gives a weird 0d numpy array
                bandList = ['u', 'g', 'r', 'i', 'z', 'y']  # TODO it better

                if patch_counter in data['patch']:
                    for band in bandList:
                        color = bands_dict['colors'][band]
                        markerstyle = bands_dict['symbols'][band]
                        bandSelector = BandSelector(vectorKey='band', bands=[band])
                        band_mask = bandSelector(data)

                        tot_mask = (patch_mask) & (band_mask)

                        ax.plot(data['depth'][tot_mask], data['pixels'][tot_mask],
                                color=color, linewidth=0, ls=None,
                                marker=markerstyle, alpha=0.75,
                                label=f'{band}')
                        ax.grid(alpha=0.5)

                    # chart formatting
                    ax.set_yscale('log')
                    ax.set_xlim(0, max_depth)
                    ax.set_ylim(5, max_pixels)
                    ax.set_xticks(np.arange(0, max_depth, 20))
                if (n != 0):
                    ax.set_yticklabels([])
                if (patch_counter not in range(10)):
                    ax.set_xticklabels([])
                ax.set_title(f"patch {patch_counter}", y=0.9)
                patch_counter += 1
                m += 1
            patch_counter -= 2*(n+1)
        fig.supxlabel('Number of input visits (n_image depth)')
        fig.supylabel('Count (pixels)')

        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
