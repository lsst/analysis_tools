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

import math
import matplotlib.pyplot as plt

from lsst.pex.config import ChoiceField, DictField, Field, FieldValidationError, ListField
from lsst.utils.plotting import make_figure, publication_plots, set_rubin_plotstyle
bands_dict = publication_plots.get_band_dicts()

from astropy.table import unique, Table
from matplotlib.ticker import SymmetricalLogLocator
from matplotlib.figure import Figure

from ...interfaces import PlotAction, Vector
from .plotUtils import addPlotInfo

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from ...interfaces import KeyedData, KeyedDataSchema


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
        # set_rubin_plotstyle()
        fig = plt.figure(dpi=300)

        patch_list = range(100)
    
        ncols = 10
        nrows = 10
    
        max_depth = max(data["depth"])
        max_pixels = max(data["pixels"])
    
        plt.figure(figsize=(4*ncols, 4*nrows))
        
        plt.subplots_adjust(hspace=0, wspace=0)
    
        patch=90 # The top left corner of a tract is patch 90
        m=0
        while patch >= 0:
            for n in range(10):
                ax = plt.subplot(nrows, ncols, m + 1)
                patch_mask = (data['patch'] == patch)
                bands = list(set(data[patch_mask]["band"]))
    
                if patch in data['patch']:
                    for band in bands:
                        color=bands_dict['colors'][band]
                        markerstyle=bands_dict['symbols'][band]
                        
                        mask = (patch_mask) & (data['band'] == band)
                        patch_band = data[mask]
                        
                        ax.plot(patch_band["depth"], patch_band["pixels"], 
                                color=color, linewidth=0, ls=None, 
                                marker=markerstyle, alpha=0.75,
                                label=f'{band}')
                        ax.grid(alpha=0.5)
        
                    # chart formatting
                    ax.set_yscale('log')
                    ax.set_xlim(0,max_depth)
                    ax.set_ylim(5,max_pixels)
                    ax.set_xticks(np.arange(0,max_depth,20))
                if (n != 0):
                    ax.set_yticklabels([])
                if (patch not in range(10)):
                    ax.set_xticklabels([])
                ax.set_title(f"patch {patch}", y=0.9)
                patch += 1
                m += 1
            patch -= 2*(n+1)
        fig.supxlabel('Number of input visits (n_image depth)')
        fig.supylabel('Count (pixels)')

        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
