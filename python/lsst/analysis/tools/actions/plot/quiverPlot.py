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

__all__ = ("QuiverPlot",)

import logging
from typing import Mapping, Optional

import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import Field
from matplotlib.figure import Figure

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from .plotUtils import addPlotInfo

_LOG = logging.getLogger(__name__)


class QuiverPlot(PlotAction):
    """Plots vectors on the detector focal plane.

    Given the posisions on the detector in x and y, the quiver
    plot draws arrows of length `length` and angle `angle`.
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis.", default="x (pixel)", optional=True)
    yAxisLabel = Field[str](doc="Label to use for the y axis.", default="y (pixel)", optional=True)
    zAxisLabel = Field[str](doc="Label to use for the arrows.", optional=True)
    qKeyLabel = Field[str](doc="Label to use for the optional quiver Key", optional=True)
    xCoordSize = Field[int]("Dimensions for X direction field to interpolate", default=4096)
    yCoordSize = Field[int]("Dimensions for Y direction field to interpolate", default=4096)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        base.append(("x", Vector))
        base.append(("y", Vector))
        base.append(("angle", Vector))
        base.append(("length", Vector))

        return base

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """NOTE currently can only check that something is not a Scalar, not
        check that the data is consistent with Vector
        """
        needed = self.getInputSchema(**kwargs)
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in input")
        for name, typ in needed:
            isScalar = issubclass((colType := type(data[name.format(**kwargs)])), Scalar)
            if isScalar and typ != Scalar:
                raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}")

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def makePlot(self, data: KeyedData, plotInfo: Optional[Mapping[str, str]] = None, **kwargs) -> Figure:

        quiverConf = {
            "pivot": "mid",
            "color": "blue",
            "width": 0.004,
        }

        dataSelector = np.isfinite(data["angle"]) & np.isfinite(data["length"])
        dataX = data["x"][dataSelector]
        dataY = data["y"][dataSelector]
        dataA = data["angle"][dataSelector]
        dataL = data["length"][dataSelector]

        U = dataL * np.cos(dataA)
        V = dataL * np.sin(dataA)

        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(111)

        q = ax.quiver(dataX, dataY, U, V, **quiverConf)
        if hasattr(self, "qKeyLabel"):
            ax.quiverkey(q, 0.9, 0.9, 1, self.qKeyLabel, labelpos="E", coordinates="figure")

        ax.set_xlim(0, self.xCoordSize)
        ax.set_ylim(0, self.yCoordSize)
        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)
        ax.set_aspect("equal", "box")

        plt.subplots_adjust(wspace=0.0, hspace=0.0, right=0.85)

        # add general plot info
        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
