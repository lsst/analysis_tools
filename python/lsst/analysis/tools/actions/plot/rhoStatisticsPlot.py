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

__all__ = ("RhoStatisticsPlot",)

from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
from matplotlib.figure import Figure

from lsst.pex.config import ConfigDictField

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from .plotUtils import addPlotInfo
from .xyPlot import XYPlot


class RhoStatisticsPlot(PlotAction):
    """Make multiple plots of rho statistics.

    Rho statistics capture the spatial correlation amongst various PSF size and
    shape residual quantities. For exact definitions, see
    :ref:`here <rho_definitions>`.
    """

    rhoPlots = ConfigDictField[str, XYPlot](
        doc="A configurable dict describing the rho statistics to plot.",
        default={},
    )

    def setDefaults(self) -> None:
        super().setDefaults()
        self.rhoPlots = {rhoName: XYPlot() for rhoName in ("rho3alt", "rho1", "rho2", "rho3", "rho4", "rho5")}

        yLabels = {
            "rho3alt": r"$\rho'_{3}(\theta) = \langle \frac{\delta T}{T}, \frac{\delta T}{T}\rangle$",
            "rho1": r"$\rho_{1}(\theta) = \langle \delta e, \delta e \rangle$",
            "rho2": r"$\rho_{2}(\theta) = \langle e, \delta e \rangle$",
            "rho3": r"$\rho_{3}(\theta) = \langle e\frac{\delta T}{T} , e\frac{\delta T}{T} \rangle$",
            "rho4": r"$\rho_{4}(\theta) = \langle \delta e, e\frac{\delta T}{T} \rangle$",
            "rho5": r"$\rho_{5}(\theta) = \langle e, e\frac{\delta T}{T} \rangle$",
        }

        for rhoId, rhoPlot in self.rhoPlots.items():
            rhoPlot.xAxisLabel = "Separation [arcmin]"
            rhoPlot.yAxisLabel = yLabels[rhoId]
            rhoPlot.xScale = "log"
            rhoPlot.yScale = "symlog"
            rhoPlot.yLinThresh = 1e-6
            rhoPlot.yLine = 0.0

        self.rhoPlots["rho3alt"].yScale = "linear"  # type: ignore

    def getInputSchema(self) -> KeyedDataSchema:
        # Docstring inherited
        base: list[tuple[str, type[Vector]]] = []
        base.append(("coord_ra", Vector))
        base.append(("coord_dec", Vector))
        base.append(("{{band}}_ixx", Vector))
        base.append(("{{band}}_iyy", Vector))
        base.append(("{{band}}_ixy", Vector))
        base.append(("{{band}}_ixxPSF", Vector))
        base.append(("{{band}}_iyyPSF", Vector))
        base.append(("{{band}}_ixyPSF", Vector))
        return base

    def getOutputNames(self) -> Iterable[str]:
        # Docstring inherited
        yield from self.rhoPlots.keys()

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure]:
        self._validateInput(data)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData) -> None:
        required = set(self.rhoPlots.keys())
        if not required.issubset(data.keys()):
            raise ValueError(f"Input data must contain {', '.join(self.rhoPlots.keys())}")

    def makePlot(
        self, data: KeyedData, plotInfo: Mapping[str, str] | None = None, **kwargs: Any
    ) -> Mapping[str, Figure]:
        r"""Make the plot(s).

        Parameters
        ----------
        data : `~pandas.core.frame.DataFrame`
            The catalog containing various rho statistics.
        plotInfo : `dict`, optional
            A dictionary of information about the data being plotted with keys:
                ``"run"``
                    The output run for the plots (`str`).
                ``"skymap"``
                    The type of skymap used for the data (`str`).
                ``"filter"``
                    The filter used for this data (`str`).
                ``"tract"``
                    The tract that the data comes from (`str`).
        **kwargs
            Additional keyword arguments to pass to the plot

        Returns
        -------
        fig_dict : `dict` [`~matplotlib.figure.Figure`]
            The resulting figures.
            The figure corresponding :math:`\rho_1(\theta)` can be accessed
            with the key `rho1` and similarly for the other rho statistics.
            :math:`\rho_3'` is accessed with the key `rho3alt`.

        Examples
        --------
        An example rho statistics plot may be seen below:

        .. image:: /_static/analysis_tools/rhoPlotExample.png

        For further details on how to generate a plot, please refer to the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """
        fig_dict: dict[str, Figure] = {}
        for rho_name in self.rhoPlots.keys():
            rho: XYPlot = self.rhoPlots[rho_name]
            subdata = {
                "x": data[rho_name].meanr,  # type: ignore
                "xerr": None,
            }
            if rho_name == "rho3alt":
                subdata["y"] = data[rho_name].xi  # type: ignore
                subdata["yerr"] = np.sqrt(data[rho_name].varxi)
            else:
                subdata["y"] = data[rho_name].xip  # type: ignore
                subdata["yerr"] = np.sqrt(data[rho_name].varxip)  # type: ignore
            fig = rho(subdata, **kwargs)
            if plotInfo is not None:
                fig_dict[rho_name] = addPlotInfo(fig, plotInfo)

        return fig_dict
