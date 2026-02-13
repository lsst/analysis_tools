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

__all__ = ("InterpolateDetectorMetricPlot",)

import logging
from collections.abc import Mapping

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from scipy.interpolate import CloughTocher2DInterpolator

from lsst.pex.config import Field, ListField

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from .plotUtils import addPlotInfo

_LOG = logging.getLogger(__name__)


class InterpolateDetectorMetricPlot(PlotAction):
    """Interpolate metrics evaluated at locations across a detector.

    The provided list of metric names and labels enables the creation of a
    multi-panel plot, with the 2D interpolation of the input metric values
    sampled on the given detector x and y coordinates.
    The interpolation evaluation grid can be controlled with the margin
    and number of grid points.
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis.", default="x (pixel)", optional=True)
    yAxisLabel = Field[str](doc="Label to use for the y axis.", default="y (pixel)", optional=True)
    zAxisLabels = ListField[str](doc="Labels to use for the z axis.", optional=True)
    metricNames = ListField[str](doc="Metrics to pull data from for interpolation", optional=False)

    xCoordSize = Field[int]("Dimensions for X direction field to interpolate", default=4096)
    yCoordSize = Field[int]("Dimensions for Y direction field to interpolate", default=4096)
    nGridPoints = Field[int]("N points in the grid for the field to interpolate", default=50)
    gridMargin = Field[int]("Grid margins for the field to interpolate", default=20)

    def getInputSchema(self) -> KeyedDataSchema:
        base = []
        base.append(("x", Vector))
        base.append(("y", Vector))
        for metricName in self.metricNames:
            base.append((metricName, Vector))

        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        return self.makePlot(data, **kwargs)

    def makePlot(self, data: KeyedData, plotInfo: Mapping[str, str] | None = None, **kwargs) -> Figure:
        """Makes a plot of a smooth interpolation of randomly
        sampled metrics in the image domain.

        Parameters
        ----------
        data : `KeyedData`
            The catalog to plot the points from, the catalog needs
            to have columns:

            * ``"x"``
                The x image coordinate of the input metric values
            * ``"y"``
                The y image coordinate of the input metric values
            * metricNames
                The column names of each image metric that needs to be
                interpolated.

        plotInfo : `dict`
            Optional. A dictionary of information about the data being plotted

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.

        Notes
        -----
        Uses the zAxisLabels config option to write the metric units and title
        for each of the used panels.
        The number of plots is determined from the number of `metricNames` in
        the config options. The colorbar of the interpolation is included for
        each panel, as well as a scatter plot showing the locations of the
        metric sampling locations.

        Examples
        --------
        An example of the plot produced from this code is here:

        .. image:: /_static/analysis_tools/interpolateDetectorPlotExample.png

        For a detailed example of how to make a plot from the command line
        please see the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """
        n_plots = len(self.metricNames)

        # boxsize has some extra space in x axis for the colorbar
        boxsize = (self.xCoordSize // (2**9), self.yCoordSize // (2**9) - 1)

        if n_plots == 1:
            fig, ax = plt.subplots(1, 1, figsize=boxsize)
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            axes = np.array([ax])
        elif n_plots <= 4:
            n_xplots, n_yplots = n_plots, 1
            fig, axes = plt.subplots(ncols=n_plots, nrows=1, figsize=(boxsize[0] * n_plots, boxsize[1]))
        else:
            # case where n_plots > 4:
            n_xplots = 4
            n_yplots = n_plots // 4 if n_plots % 4 == 0 else n_plots // 4 + 1
            fig, axes = plt.subplots(
                ncols=n_xplots,
                nrows=n_yplots,
                figsize=(boxsize[0] * n_xplots, boxsize[1] * n_yplots),
                sharex=True,
                sharey=True,
            )

        ytox_ratio = self.yCoordSize // self.xCoordSize
        X = np.linspace(-self.gridMargin, self.xCoordSize + self.gridMargin, self.nGridPoints)
        Y = np.linspace(-self.gridMargin, self.yCoordSize + self.gridMargin, self.nGridPoints * ytox_ratio)
        meshgridX, meshgridY = np.meshgrid(X, Y)  # 2D grid for interpolation

        if self.zAxisLabels is None:
            zAxisLabels = ["px_frac" for metricName in self.metricNames]
        else:
            zAxisLabels = self.zAxisLabels

        for ax, metricName, zlabel in zip(axes.flatten(), self.metricNames, zAxisLabels):
            dataSelector = np.isfinite(data[metricName])
            if np.count_nonzero(dataSelector) < 4:
                # Need at least four valid points for the interpolation.
                ax.set_aspect("equal", "box")
                ax.set_title(f"{metricName}[{zlabel}]")

                continue

            dataX = data["x"][dataSelector]
            dataY = data["y"][dataSelector]
            dataZ = data[metricName][dataSelector]

            interp = CloughTocher2DInterpolator(list(zip(dataX, dataY)), dataZ)
            Z = interp(meshgridX, meshgridY)

            pc = ax.pcolormesh(X, Y, Z, shading="auto")
            ax.scatter(dataX, dataY, s=10, facecolor="silver", edgecolor="black")
            _ = fig.colorbar(pc, shrink=0.7, location="right", fraction=0.07)
            ax.set_aspect("equal", "box")
            ax.set_title(f"{metricName}[{zlabel}]")

        for iax in range(axes.size - n_plots):
            axes.flatten()[-(iax + 1)].remove()
        # setting x axis labels for all
        # setting only y axis labels for plots on the left of the panel
        if n_yplots > 1:
            for ax in axes[:, 0]:
                ax.set_ylabel("Y")
            for ax in axes[-1, :]:
                ax.set_xlabel("X")
        else:
            axes[0].set_ylabel("Y")
            for ax in axes:
                ax.set_xlabel("X")

        plt.tight_layout()

        # add general plot info
        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
