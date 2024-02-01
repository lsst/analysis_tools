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

__all__ = ("GridPlot", "GridPanelConfig")

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from lsst.pex.config import Config, ConfigDictField, DictField, Field, ListField
from lsst.pex.config.configurableActions import ConfigurableActionField
from matplotlib.gridspec import GridSpec

from ...interfaces import PlotAction, PlotElement

if TYPE_CHECKING:
    from lsst.analysis.tools.interfaces import KeyedData, PlotResultType


class GridPanelConfig(Config):
    plotElement = ConfigurableActionField[PlotElement](
        doc="Plot element.",
    )
    title = DictField[str, str](
        doc="String arguments passed into ax.set_title() defining the plot element title.",
    )
    titleY = Field[float](
        doc="Y position of plot element title.",
        default=None,
    )


class GridPlot(PlotAction):
    """Plot a series of plot elements onto a regularly spaced grid."""

    plotElements = ConfigDictField(
        doc="Plot elements.",
        keytype=int,
        itemtype=GridPanelConfig,
    )
    numRows = Field[int](
        doc="Number of rows.",
        default=1,
    )
    numCols = Field[int](
        doc="Number of columns.",
        default=1,
    )
    xDataKeys = DictField[int, str](
        doc="Dependent data definitions. The key of this dict is the panel ID. The values are keys of data "
        "to plot (comma-separated for multiple) where each key may be a subset of a full key.",
        default={},
    )
    valsGroupBy = DictField[int, str](
        doc="Independent data definitions. The key of this dict is the panel ID. The values are keys of data "
        "to plot (comma-separated for multiple) where each key may be a subset of a full key.",
    )
    figsize = ListField[float](
        doc="Figure size.",
        default=[8, 8],
    )
    dpi = Field[float](
        doc="Dots per inch.",
        default=150,
    )
    suptitle = DictField[str, str](
        doc="String arguments passed into fig.suptitle() defining the figure title.",
        optional=True,
    )

    def __call__(self, data: KeyedData, **kwargs) -> PlotResultType:
        """Plot data."""
        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        if self.suptitle is not None:
            fig.suptitle(**self.suptitle)
        gs = GridSpec(self.numRows, self.numCols, figure=fig)

        for row in range(self.numRows):
            for col in range(self.numCols):
                index = row * self.numCols + col
                if index not in self.valsGroupBy.keys():
                    continue
                ax = fig.add_subplot(gs[row, col])

                xList = x.split(",") if (x := self.xDataKeys.get(index)) else None
                valList = self.valsGroupBy[index].split(",")

                for i, val in enumerate(valList):
                    for key in data:
                        newData = {}
                        if val not in key:
                            continue
                        namedKey = self.plotElements[index].plotElement.valsKey
                        newData[namedKey] = data[key]
                        if xList is not None:
                            namedKey = self.plotElements[index].plotElement.xKey
                            newData[namedKey] = data[xList[i]]

                        _ = self.plotElements[index].plotElement(data=newData, ax=ax, **kwargs)

                if self.plotElements[index].title is not None:
                    ax.set_title(**self.plotElements[index].title, y=self.plotElements[index].titleY)

        plt.tight_layout()
        fig.show()
        return fig

    def validate(self):
        """Validate configuration."""
        super().validate()
        if self.xDataKeys and len(self.xDataKeys) != self.numRows * self.numCols:
            raise RuntimeError("Number of xDataKeys keys must match number of rows * columns.")
        if len(self.valsGroupBy) != self.numRows * self.numCols:
            raise RuntimeError("Number of valsGroupBy keys must match number of rows * columns.")
