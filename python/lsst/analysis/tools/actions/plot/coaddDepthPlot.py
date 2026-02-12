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

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from lsst.skymap.tractInfo import ExplicitTractInfo
from lsst.utils.plotting import publication_plots, set_rubin_plotstyle

from ...interfaces import PlotAction, Vector
from .plotUtils import addPlotInfo

if TYPE_CHECKING:
    from ...interfaces import KeyedData, KeyedDataSchema

bands_dict = publication_plots.get_band_dicts()


class CoaddDepthPlot(PlotAction):
    """Make a plot of pixels per coadd depth."""

    def setDefaults(self):
        super().setDefaults()

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector]]] = []
        base.append(("patch", Vector))
        base.append(("band", Vector))
        base.append(("depth", Vector))
        base.append(("pixels", Vector))
        return base

    def __call__(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        **kwargs,
    ) -> Figure:
        self._validateInput(data, tractInfo)

        if not isinstance(tractInfo, ExplicitTractInfo):
            raise TypeError(f"Input `tractInfo` type must be {ExplicitTractInfo} not {type(tractInfo)}.")

        return self.makePlot(data, tractInfo, **kwargs)

    def _validateInput(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
    ) -> None:
        needed = set(k[0] for k in self.getInputSchema())
        if not needed.issubset(data.keys()):
            raise ValueError(f"Input data does not contain all required keys: {self.getInputSchema()}")

    def makePlot(
        self,
        data: KeyedData,
        tractInfo: ExplicitTractInfo,
        plotInfo: Mapping[str, str] | None = None,
        **kwargs: Any,
    ) -> Figure:
        """Make the plot.

        Parameters
        ----------
        `KeyedData`
            The catalog to plot the points from.
        tractInfo : `~lsst.skymap.tractInfo.ExplicitTractInfo`
            The tract info object.
        plotInfo : `dict`
            A dictionary of the plot information.

        Returns
        -------
        fig : `~matplotlib.figure.Figure`
            The resulting figure.

        Examples
        --------
        An example coaddDepthPlot may be seen below:

        .. image:: /_static/analysis_tools/coaddDepthPlotExample.png
        """
        set_rubin_plotstyle()
        fig = plt.figure(dpi=300, figsize=(20, 20))

        # Get number of vertical and horizontal patches in the tract.
        num_patches_x, num_patches_y = tractInfo.getNumPatches()

        max_depth = max(data["depth"])
        max_pixels = max(data["pixels"])

        plt.subplots_adjust(hspace=0, wspace=0)

        patch_counter = (num_patches_x * num_patches_y) - num_patches_x  # The top left corner of a tract
        m = 0  # subplot index
        while patch_counter >= 0:
            for n in range(num_patches_x):  # column index
                ax = plt.subplot(num_patches_x, num_patches_y, m + 1)
                patch_mask = data["patch"] == patch_counter

                if patch_counter in data["patch"]:
                    uniqueBands = set(data["band"][patch_mask])
                    for band in uniqueBands:
                        color = bands_dict["colors"][band]
                        markerstyle = bands_dict["symbols"][band]
                        band_mask = data["band"] == band

                        tot_mask = (patch_mask) & (band_mask)

                        ax.plot(
                            data["depth"][tot_mask],
                            data["pixels"][tot_mask],
                            color=color,
                            linewidth=0,
                            ls=None,
                            marker=markerstyle,
                            ms=4,
                            alpha=0.75,
                            label=f"{band}",
                        )
                        ax.grid(alpha=0.5)

                # Chart formatting
                # Need a solution for ax.set_xscale when max_depth is high,
                # but not all patches/bands have quality coverage.
                ax.set_yscale("log")
                ax.set_xlim(0, max_depth + 5)
                ax.set_ylim(5, max_pixels)
                # Can we somehow generalize ax.set_xticks?
                # ax.set_xticks(np.arange(0, max_depth, 20))
                ax.tick_params(axis="both", which="minor")
                ax.tick_params(axis="both", which="both", top=False, right=False)

                # Only label axes of the farmost left and bottom row of plots
                if n != 0:
                    ax.set_yticklabels([])
                    ax.tick_params(axis="y", which="both", length=0)
                if patch_counter not in range(num_patches_x):
                    ax.set_xticklabels([])
                    ax.tick_params(axis="x", which="both", length=0)

                ax.set_title(f"patch {patch_counter}", y=0.85)
                patch_counter += 1
                m += 1
            patch_counter -= 2 * (n + 1)
        fig.supxlabel("Number of input visits (n_image depth)", y=0.075)
        fig.supylabel("Count (pixels)", x=0.075)
        legend_elements = [
            Line2D(
                [0], [0], color=bands_dict["colors"]["u"], lw=0, marker=bands_dict["symbols"]["u"], label="u"
            ),
            Line2D(
                [0], [0], color=bands_dict["colors"]["g"], lw=0, marker=bands_dict["symbols"]["g"], label="g"
            ),
            Line2D(
                [0], [0], color=bands_dict["colors"]["r"], lw=0, marker=bands_dict["symbols"]["r"], label="r"
            ),
            Line2D(
                [0], [0], color=bands_dict["colors"]["i"], lw=0, marker=bands_dict["symbols"]["i"], label="i"
            ),
            Line2D(
                [0], [0], color=bands_dict["colors"]["z"], lw=0, marker=bands_dict["symbols"]["z"], label="z"
            ),
            Line2D(
                [0], [0], color=bands_dict["colors"]["y"], lw=0, marker=bands_dict["symbols"]["y"], label="y"
            ),
        ]
        plt.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(1.05, 10))

        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
