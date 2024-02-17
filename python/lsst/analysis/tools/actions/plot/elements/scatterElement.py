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

__all__ = ("ScatterElement",)

from typing import TYPE_CHECKING

from lsst.pex.config import Field

from ....interfaces import PlotElement

if TYPE_CHECKING:
    from lsst.analysis.tools.interfaces import KeyedData
    from matplotlib.axes import Axes


class ScatterElement(PlotElement):
    """Configuration options for scatter plot plot elements.

    Attributes
    ----------
    xKey : `~lsst.pex.config.Field`
        X-axis data vector key. If None, the index of the values are used.
    valsKey : `~lsst.pex.config.Field`
        Y-axis data vector key. Plotted against the index if no xKey is given.
    color : `~lsst.pex.config.Field`
        Color.
    marker : `~lsst.pex.config.Field`
        Point marker.
    linestyle : `~lsst.pex.config.Field`
        Linestyle.
    linewidth : `~lsst.pex.config.Field`
        Linewidth.
    markersize : `~lsst.pex.config.Field`
        Markersize.
    """

    xKey = Field[str](
        doc="X-axis data vector key. If None, the index of the values are used.",
        optional=True,
    )
    valsKey = Field[str](
        doc="Y-axis data vector key. Plotted against the index if no xKey is given.",
        default="value",
    )
    color = Field[str](
        doc="Color",
        optional=True,
    )
    marker = Field[str](
        doc="Point marker",
        optional=True,
    )
    linestyle = Field[str](
        doc="Linestyle",
        optional=True,
    )
    linewidth = Field[float](
        doc="Linewidth",
        optional=True,
    )
    markersize = Field[float](
        doc="Markersize",
        optional=True,
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
        self._validateInputs(data)
        ax.plot(
            data[self.xKey] if self.xKey is not None else range(len(data[self.valsKey])),  # type: ignore
            data[self.valsKey],  # type: ignore
            color=self.color if self.color is not None else None,
            marker=self.marker if self.marker is not None else None,
            linestyle=self.linestyle if self.linestyle is not None else None,
            linewidth=self.linewidth if self.linewidth is not None else None,
            markersize=self.markersize if self.markersize is not None else None,
        )

        return data

    def _validateInputs(self, data: KeyedData) -> None:
        """Validate inputs.

        Parameters
        ----------
        data : `~lsst.analysis.tools.interfaces.KeyedData`
            Data to plot.
        """
        if self.xKey is not None and len(data[self.xKey]) != len(data[self.valsKey]):  # type: ignore
            raise ValueError("X and Y vector inputs must be the same length.")
