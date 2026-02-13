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

__all__ = ("DiaSkyPanel", "DiaSkyPlot")

from collections.abc import Mapping

import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from lsst.pex.config import ConfigDictField, Field, ListField

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from .plotUtils import PanelConfig


class DiaSkyPanel(PanelConfig):
    """Configuration options for DiaSkyPlot panels."""

    xlabel = Field[str](
        doc="Panel x-axis label.",
        default="RA (deg)",
    )
    ylabel = Field[str](
        doc="Panel y-axis label.",
        default="Dec (deg)",
    )
    invertXAxis = Field[bool](
        doc="Invert x-axis?",
        default=True,
    )
    size = Field[float](
        doc="Point size",
        default=20,
    )
    alpha = Field[float](
        doc="Point transparency",
        default=0.5,
    )
    # Eventually we might retrieve data from more columns to make the plot
    # prettier/more information rich
    ras = ListField[str](
        doc="Names of RA columns",
        optional=False,
    )
    decs = ListField[str](
        doc="Names of Dec columns",
        optional=False,
    )
    colorList = Field[str](
        doc="Colors for the points",
        optional=True,
    )
    legendLabels = ListField[str](
        doc="Labels for the legend",
        optional=True,
    )


class DiaSkyPlot(PlotAction):
    """Generic pseudo base class for plotting DiaSources
    (or DiaObjects) on the sky.
    """

    panels = ConfigDictField(
        doc="A configurable dict describing the panels to be plotted (both data columns and layouts).",
        keytype=str,
        itemtype=DiaSkyPanel,
        default={},
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        """Defines the schema this plot action expects (the keys it looks
        for and what type they should be). In other words, verifies that
        the input data has the columns we are expecting with the right dtypes.
        """
        for ra in self.panels.ras.values():
            yield (ra, Vector)
        for dec in self.panels.decs.values():
            yield (dec, Vector)

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        return self.makePlot(data, **kwargs)

    def makePlot(self, data: KeyedData, **kwargs) -> Figure:
        """Make an N-panel plot with locations of DiaSources or
        DiaObjects displayed in each panel.

        Parameters
        ----------
        data : `lsst.analysis.tools.interfaces.KeyedData`

        Returns
        -------
        fig : `matplotlib.figure.Figure`
        """
        if "figsize" in kwargs:
            figsize = kwargs.pop("figsize", "")
            fig = plt.figure(figsize=figsize, dpi=300)
        else:
            fig = plt.figure(figsize=(12, 9), dpi=300)
        axs = self._makeAxes(fig)
        for panel, ax in zip(self.panels.values(), axs):
            self._makePanel(data, panel, ax, **kwargs)
        plt.draw()
        return fig

    def _makeAxes(self, fig):
        """Determine axes layout for main figure.

        Use matplotlib's subplot2grid to determine the panel geometry,
        which calls gridspec.

        Parameters
        ----------
        fig : `matplotlib.figure.Figure`

        Returns
        -------
        axs : `list` containing one or more matplotlib axes, one for each panel
        """
        axs = []
        for count, panel in enumerate(self.panels.values()):
            subplot2gridShape = (panel.subplot2gridShapeRow, panel.subplot2gridShapeColumn)
            subplot2gridLoc = (panel.subplot2gridLocRow, panel.subplot2gridLocColumn)
            axs.append(
                plt.subplot2grid(
                    subplot2gridShape,
                    subplot2gridLoc,
                    rowspan=panel.subplot2gridRowspan,
                    colspan=panel.subplot2gridColspan,
                )
            )
        return axs

    def _makePanel(self, data, panel, ax, **kwargs):
        """Plot a single panel.

        Parameters
        ----------
        data : `lsst.analysis.tools.interfaces.KeyedData`
        panel : `DiaSkyPanel`
        ax : matplotlib axis
        color : `str`
        """
        artists = []  # Placeholder for each series being plotted
        for idx, (ra, dec) in enumerate(zip(panel.ras, panel.decs)):  # loop over column names (dict keys)
            if panel.colorList:
                color = panel.colorList[idx]
                artist = ax.scatter(
                    data[ra], data[dec], s=panel.size, alpha=panel.alpha, marker=".", linewidths=0, c=color
                )
            else:  # Use matplotlib default colors
                artist = ax.scatter(
                    data[ra], data[dec], s=panel.size, alpha=panel.alpha, marker=".", linewidths=0
                )
            artists.append(artist)
            # TODO DM-42768: implement lists of sizes, alphas, etc.
            # and add better support for multi-panel plots.

        ax.set_xlabel(panel.xlabel)
        ax.set_ylabel(panel.ylabel)
        if panel.legendLabels:
            ax.legend(artists, panel.legendLabels)
        if panel.invertXAxis:
            ax.invert_xaxis()
        if panel.topSpinesVisible:
            ax.spines["top"].set_visible(True)
        else:
            ax.spines["top"].set_visible(False)
        if not panel.bottomSpinesVisible:  # default is True
            ax.spines["bottom"].set_visible(False)
        if not panel.leftSpinesVisible:
            # Default is True; if False, also put ticks and labels on the right
            ax.spines["left"].set_visible(False)
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")
        if not panel.rightSpinesVisible:  # default is True
            ax.spines["right"].set_visible(False)
