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

__all__ = ("PlotInfoElement",)

from typing import TYPE_CHECKING, Any, Mapping

import matplotlib.transforms as transforms

from ....interfaces import PlotElement
from ..plotUtils import parsePlotInfo

if TYPE_CHECKING:
    from matplotlib.axes import Axes


class PlotInfoElement(PlotElement):
    """Adds information to a plot, such as run collection name,
    datasets used, dataset type, and tract/visit IDs.

    Notes
    -----
    Converts a plotInfo dict to a string and creates a textbox to display it.

    """

    def __call__(self, ax: Axes, figureInfo: Mapping[str, Any], **kwargs) -> Axes:
        """Draw a text box that displays useful information about the plot.

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            Axes to plot on.
        figInfo: `dict`
            A dictionary containing information about the figure on which the
            element is being placed.  It should have the following keys:

            ``"figsize"``
                The width and height of the figure in inches
                (`tuple` of two `float`s)
            ``"dpi"``
                The dpi of the figure (`float`)

        Notes
        -----
        plotInfoElement is added to the figure in-place, so no parameter is
        returned.
        """
        plotInfo = kwargs["plotInfo"]
        plotInfoText = parsePlotInfo(plotInfo)

        # Ensures consistent spacing between the plotName and plotInfoText:
        dpi_scale_trans = transforms.Affine2D().scale(figureInfo["dpi"])
        trans = ax.transAxes + transforms.ScaledTranslation(0.0, -0.28, dpi_scale_trans)

        ax.set_axis_off()
        ax.text(0, 1.0, plotInfo["plotName"], fontsize=7, ha="left", va="bottom")
        ax.text(0, 1.0, plotInfoText, fontsize=6, alpha=0.6, ha="left", transform=trans)
