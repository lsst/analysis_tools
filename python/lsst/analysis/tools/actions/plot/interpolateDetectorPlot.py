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
from typing import Mapping, Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy.interpolate import CloughTocher2DInterpolator

from lsst.pex.config import Field

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from .plotUtils import addPlotInfo

_LOG = logging.getLogger(__name__)


class InterpolateDetectorMetricPlot(PlotAction):
    """Returns a fixed array size with interpolation of metric in the detector
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis.", default="x (pixel)", optional=True)
    yAxisLabel = Field[str](doc="Label to use for the y axis.", default="y (pixel)", optional=True)
    zAxisLabel = Field[str](doc="Label to use for the z axis.", optional=True)

    xCoordSize = Field[int]("Dimensions for X direction field to interpolate", default=4000)
    yCoordSize = Field[int]("Dimensions for Y direction field to interpolate", default=4072)
    nGridPoints = Field[int]("N points in the grid for the field to interpolate", default=40)
    gridMargin = Field[int]("Grid margins for the field to interpolate", default=20)

    def getInputSchema(self) -> KeyedDataSchema:
        base = []

        base.append(("x", Vector))
        base.append(("y", Vector))
        base.append(("metricValues", Vector))

        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        # self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def makePlot(
            self,
            data: KeyedData,
            plotInfo: Optional[Mapping[str, str]] = None,
            **kwargs
    ) -> Figure:

        X = np.linspace(-self.gridMargin, self.xCoordSize+self.gridMargin, self.nGridPoints)
        Y = np.linspace(-self.gridMargin, self.yCoordSize+self.gridMargin, self.nGridPoints)
        meshgridX, meshgridY = np.meshgrid(X, Y)  # 2D grid for interpolation

        interp = CloughTocher2DInterpolator(list(zip(data["x"], data["y"])), data["metricValues"])
        Z = interp(meshgridX, meshgridY)

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        pc = ax.pcolormesh(X, Y, Z, shading='auto')
        ax.scatter(data['x'], data['y'], s=5, c="black")
        cbar = fig.colorbar(pc)
        cbar.set_label(self.zAxisLabel, rotation=270)
        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)
        ax.set_aspect('equal', 'box')
        
        # add general plot info
        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
