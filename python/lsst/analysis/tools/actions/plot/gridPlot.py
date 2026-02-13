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
import numpy as np
from matplotlib.gridspec import GridSpec

from lsst.pex.config import Config, ConfigDictField, DictField, Field, ListField
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import PlotAction, PlotElement

if TYPE_CHECKING:
    from lsst.analysis.tools.interfaces import KeyedData, PlotResultType


class GridPanelConfig(Config):
    plotElement = ConfigurableActionField[PlotElement](
        doc="Plot element.",
    )
    title = DictField[str, str](
        doc="String arguments passed into ax.set_title() defining the plot element title.",
        optional=True,
    )
    titleY = Field[float](
        doc="Y position of plot element title.",
        optional=True,
    )


class GridPlot(PlotAction):
    """Plot a series of plot elements onto a regularly spaced grid."""

    panels = ConfigDictField(
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
    width_ratios = ListField[float](
        doc="Width ratios",
        optional=True,
    )
    height_ratios = ListField[float](
        doc="Height ratios",
        optional=True,
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
    xAxisLabel = Field[str](
        doc="String argument passed into fig.supxlabel() defining the figure x label.",
        optional=True,
    )
    yAxisLabel = Field[str](
        doc="String argument passed into fig.supylabel() defining the figure y label.",
        optional=True,
    )

    def __call__(self, data: KeyedData, **kwargs) -> PlotResultType:
        """Plot data."""
        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        figureInfo = {"figsize": fig.get_size_inches(), "dpi": fig.get_dpi()}

        if self.height_ratios is None:
            height_ratios = np.ones(self.numRows) / self.numRows
        else:
            height_ratios = self.height_ratios / np.sum(self.height_ratios)

        if self.width_ratios is None:
            width_ratios = np.ones(self.numCols) / self.numCols
        else:
            width_ratios = self.width_ratios / np.sum(self.width_ratios)

        if self.suptitle is not None:
            fig.suptitle(**self.suptitle)
        if self.xAxisLabel is not None:
            fig.supxlabel(self.xAxisLabel)
        if self.yAxisLabel is not None:
            fig.supylabel(self.yAxisLabel)

        # TODO: See DM-44283:Add subplot_mosaic functionality to plotElements
        gs = GridSpec(
            self.numRows,
            self.numCols,
            figure=fig,
            height_ratios=height_ratios,
            width_ratios=width_ratios,
        )

        # Iterate over all of the plots we'll make:
        for row in range(self.numRows):
            for col in range(self.numCols):
                # This sequential index is used to identify what data
                # to plot.  The `valsGroupBy` dict should have this
                # index as a key, with the values matching the subset
                # of rows that have that value in the column specified
                # by the `panelKey`.
                index = row * self.numCols + col
                if index not in self.valsGroupBy.keys():
                    continue
                ax = fig.add_subplot(gs[row, col])

                # These lists hold the columns that will be plotted,
                # comma separated to allow multiple series to be
                # plotted on the same panel.  If `xDataKeys` does not
                # contain this panel's index, then the vector index
                # will be used for the x-coordinate.
                xList = x.split(",") if (x := self.xDataKeys.get(index)) else None
                valList = self.valsGroupBy[index].split(",")

                # Iterate over the series to plot in this panel:
                for i, val in enumerate(valList):
                    for key in data:
                        newData = {}
                        if val not in key:
                            # Skip columns in data that do not match
                            # our series identifier.
                            continue
                        if xList is not None:
                            # Store the x-coordinate data to be
                            # plotted in the temporary column name
                            # indicated by the `xDataKeys` dict above.
                            namedKey = self.panels[index].plotElement.xKey
                            newData[namedKey] = data[xList[i]]
                            if key in xList:
                                # if this key is in the xList, we need
                                # to not plot it.
                                continue

                        # If provided, store the y-coordinate data to be
                        # plotted in the temporary column name indicated
                        # by the `valsGroupBy` dict above. Not all elements
                        # need y-coordinate data, such as plotInfoElement.
                        if hasattr(self.panels[index].plotElement, "valsKey"):
                            namedKey = self.panels[index].plotElement.valsKey
                            newData[namedKey] = data[key]

                        # Actually make the plot.
                        _ = self.panels[index].plotElement(
                            data=newData, ax=ax, figureInfo=figureInfo, **kwargs
                        )

                if self.panels[index].title is not None:
                    ax.set_title(**self.panels[index].title, y=self.panels[index].titleY)

        plt.tight_layout()
        return fig

    def validate(self):
        """Validate configuration."""
        super().validate()
        if self.xDataKeys and len(self.xDataKeys) != self.numRows * self.numCols:
            raise RuntimeError("Number of xDataKeys keys must match number of rows * columns.")
        if len(self.valsGroupBy) != self.numRows * self.numCols:
            raise RuntimeError("Number of valsGroupBy keys must match number of rows * columns.")
        if self.width_ratios and len(self.width_ratios) != self.numCols:
            raise RuntimeError("Number of supplied width ratios must match number of columns.")
        if self.height_ratios and len(self.height_ratios) != self.numRows:
            raise RuntimeError("Number of supplied height ratios must match number of rows.")
