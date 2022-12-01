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

from typing import TYPE_CHECKING, Any, Mapping

__all__ = ("XYPlot",)

import matplotlib.pyplot as plt
from lsst.pex.config import ChoiceField, DictField, Field, FieldValidationError

from ...interfaces import PlotAction, Vector
from .plotUtils import addPlotInfo

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from ...interfaces import KeyedData, KeyedDataSchema


class XYPlot(PlotAction):
    boolKwargs = DictField[str, bool](
        doc="Keyword arguments to ax.errorbar that take boolean values",
        default={},
        optional=True,
    )

    numKwargs = DictField[str, float](
        doc="Keyword arguments to ax.errorbar that take numerical (float or int) values",
        default={},
        optional=True,
    )

    strKwargs = DictField[str, str](
        doc="Keyword arguments to ax.errorbar that take string values",
        default={},
        optional=True,
    )

    xAxisLabel = Field[str](
        doc="The label to use for the x-axis.",
        default="x",
    )

    yAxisLabel = Field[str](
        doc="The label to use for the y-axis.",
        default="y",
    )

    xScale = ChoiceField[str](
        doc="The scale to use for the x-axis.",
        default="linear",
        allowed={scale: scale for scale in ("linear", "log", "symlog")},
    )

    yScale = ChoiceField[str](
        doc="The scale to use for the y-axis.",
        default="linear",
        allowed={scale: scale for scale in ("linear", "log", "symlog")},
    )

    xLinThresh = Field[float](
        doc=(
            "The value around zero where the scale becomes linear in x-axis "
            "when symlog is set as the scale. Sets the `linthresh` parameter "
            "of `~matplotlib.axes.set_xscale`."
        ),
        default=1e-6,
        optional=True,
    )

    yLinThresh = Field[float](
        doc=(
            "The value around zero where the scale becomes linear in y-axis "
            "when symlog is set as the scale. Sets the `linthresh` parameter "
            "of `~matplotlib.axes.set_yscale`."
        ),
        default=1e-6,
        optional=True,
    )

    def setDefaults(self):
        super().setDefaults()
        self.strKwargs = {"fmt": "o"}

    def validate(self):
        if (len(set(self.boolKwargs.keys()).intersection(self.numKwargs.keys())) > 0) or (
            len(set(self.boolKwargs.keys()).intersection(self.strKwargs.keys())) > 0
        ):
            raise FieldValidationError(self.boolKwargs, self, "Keywords have been repeated")

        super().validate()

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector]]] = []
        base.append(("x", Vector))
        base.append(("y", Vector))
        base.append(("xerr", Vector))
        base.append(("yerr", Vector))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Figure:
        self._validateInput(data)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData) -> None:
        needed = set(k[0] for k in self.getInputSchema())
        if not needed.issubset(data.keys()):
            raise ValueError(f"Input data does not contain all required keys: {self.getInputSchema()}")

    def makePlot(self, data: KeyedData, plotInfo: Mapping[str, str] | None = None, **kwargs: Any) -> Figure:
        """Make the plot.

        Parameters
        ----------
        data : `~pandas.core.frame.DataFrame`
            The catalog containing various rho statistics.
        **kwargs
            Additional keyword arguments to pass to the plot
        """
        # Allow for multiple curves to lie on the same plot.
        fig = kwargs.get("fig", None)
        if fig is None:
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(111)
        else:
            ax = fig.gca()

        ax.errorbar(
            data["x"],
            data["y"],
            xerr=data["xerr"],
            yerr=data["yerr"],
            **self.boolKwargs,  # type: ignore
            **self.numKwargs,  # type: ignore
            **self.strKwargs,  # type: ignore
        )
        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)
        if self.xScale == "symlog":
            ax.set_xscale("symlog", linthresh=self.xLinThresh)
            ax.fill_betweenx(y=data["y"], x1=-self.xLinThresh, x2=self.xLinThresh, color="gray", alpha=0.2)
        else:
            ax.set_xscale(self.xScale)  # type: ignore
        if self.yScale == "symlog":
            ax.set_yscale("symlog", linthresh=self.yLinThresh)
            ax.fill_between(x=data["x"], y1=-self.yLinThresh, y2=self.yLinThresh, color="gray", alpha=0.2)
        else:
            ax.set_yscale(self.yScale)  # type: ignore

        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
