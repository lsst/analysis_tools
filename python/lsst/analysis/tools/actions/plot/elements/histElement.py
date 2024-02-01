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

__all__ = ("HistElement",)

from typing import TYPE_CHECKING

from lsst.pex.config import Field

from ....interfaces import PlotElement

if TYPE_CHECKING:
    from lsst.analysis.tools.interfaces import KeyedData
    from matplotlib.axes import Axes


class HistElement(PlotElement):
    """Configuration options for ax.hist elements.

    Attributes
    ----------
    valsKey : `~lsst.pex.config.Field`
        Y-axis data vector key.
    bins : `~lsst.pex.config.Field`
        Number of x axis bins within plot x-range.
    color : `~lsst.pex.config.Field`
        Color.
    density : `~lsst.pex.config.Field`
        If True, draw and return a probability density.
    cumulative : `~lsst.pex.config.Field`
        If True, draw a cumulative histogram of the data. If ``density`` is
        also True, the histogram is normalized such that the last bin equals 1.
    """

    valsKey = Field[str](
        doc="Y-axis data vector key.",
        default="value",
    )
    bins = Field[int](
        doc="Number of x axis bins within plot x-range.",
        default=50,
    )
    color = Field[str](
        doc="Color",
        optional=True,
    )
    density = Field[bool](
        doc="If True, draw and return a probability density.",
        default=False,
    )
    cumulative = Field[bool](
        doc="If True, draw a cumulative histogram of the data. If ``density`` "
        "is also True, the histogram is normalized such that the last bin "
        "equals 1.",
        default=False,
    )

    def __call__(self, data: KeyedData, ax: Axes, **kwargs) -> KeyedData:
        """Plot y versus x as lines and/or markers.

        Parameters
        ----------
        data : `~lsst.analysis.tools.interfaces.KeyedData`
            Keyed data containing the data to plot.
        ax : `~matplotlib.axes.Axes`
            Axes to plot on.

        Returns
        -------
        data : `~lsst.analysis.tools.interfaces.KeyedData`
            Data used for plotting.
        """
        ax.hist(
            x=data[self.valsKey],  # type: ignore
            bins=self.bins,
            color=self.color,
            density=self.density,
            cumulative=self.cumulative,
        )

        return data
