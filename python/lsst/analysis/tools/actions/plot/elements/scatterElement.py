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

from lsst.pex.config import Field, ListField

from ....interfaces import PlotElement

if TYPE_CHECKING:
    from lsst.analysis.tools.interfaces import KeyedData
    from matplotlib.axes import Axes


class ScatterElement(PlotElement):
    """Configuration options for ScatterPlot elements.

    Attributes
    ----------
    x : `~lsst.pex.config.Field`
        X-axis data vector.
    y : `~lsst.pex.config.Field`
        Y-axis data vector.
    color : `~lsst.pex.config.Field`
        Point color.
    size : `~lsst.pex.config.Field`
        Point size.
    alpha : `~lsst.pex.config.Field`
        Point transparency.
    symbol : `~lsst.pex.config.Field`
        Point symbol.
    """

    xKey = Field[str](
        doc="X-axis data vector.",
        optional=True,
    )
    valsKey = Field[str](
        doc="Y-axis data vector.",
        default="value",
    )
    color = Field[str](
        doc="Point color",
        optional=True,
    )
    size = Field[float](
        doc="Point size",
        default=5,
    )
    alpha = Field[float](
        doc="Point transparency",
        default=1,
    )
    symbol = Field[str](
        doc="Point symbol",
        default="o",
    )

    def __call__(self, data: KeyedData, ax: Axes, **kwargs) -> KeyedData:
        """Plot data as a scatter plot.

        Parameters
        ----------
        data : `~lsst.analysis.tools.interfaces.KeyedData`
            Data to plot.
        ax : `~matplotlib.axes.Axes`
            Axes to plot on.

        Returns
        -------
        data : `~lsst.analysis.tools.interfaces.KeyedData`
            Metadata generated during plotting.
        """
        self._validateInputs(data)
        ax.plot(
            data[self.xKey] if self.xKey is not None else range(len(data[self.valsKey])),  # type: ignore
            data[self.valsKey],  # type: ignore
            # c=self.color,  # type: ignore
            # s=self.size,  # type: ignore
            # alpha=self.alpha,  # type: ignore
            # marker=self.symbol,  # type: ignore
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
