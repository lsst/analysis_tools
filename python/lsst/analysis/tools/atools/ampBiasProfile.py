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

from lsst.analysis.tools.interfaces._interfaces import KeyedDataSchema

__all__ = ("AmpBiasProfileTool",)

from lsst.pex.config import Field

from ..actions.plot.elements.scatterElement import ScatterElement
from ..actions.plot.gridPlot import GridPlot, PlotElementConfig
from ..interfaces import AnalysisTool, Vector


class AmpBiasProfileTool(AnalysisTool):
    parameterizedBand: bool = False

    panelKey = Field[str](
        doc="Panel selector.",
    )
    dataKey = Field[str](
        doc="Data selector.",
    )
    quantityKey = Field[str](
        doc="Quantity selector.",
    )

    def _call_single(self, data: KeyedData, **kwargs) -> PlotResultType:
        # Repack data here.
        newData = {}
        for i in range(len(data[self.panelKey])):
            newData[f"{data[self.panelKey][i]}_{data[self.dataKey][i]}_{self.quantityKey}"] = data[
                self.quantityKey
            ][i]
        return super()._call_single(newData, **kwargs)

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.panelKey, Vector),
            (self.dataKey, Vector),
            (self.quantityKey, Vector),
        )

    def setDefaults(self):
        super().setDefaults()
        self.panelKey = "amplifier"
        self.dataKey = "mjd"
        self.quantityKey = "biasSerialProfile"

        # # Make these input data columns available for subsequent use
        # self.prep.keysToLoad = ["amplifier", "mjd", "biasSerialProfile"]

        self.produce.plot = GridPlot()

        self.produce.plot.plotElements = {}

        self.produce.plot.numRows = 2
        self.produce.plot.numCols = 2

        self.produce.plot.valsGroupBy = {
            0: "C00",
            1: "C01",
            2: "C02",
            3: "C03",
        }
        # self.produce.plot.valsGroupBy[0] = "C00"
        # self.produce.plot.valsGroupBy[1] = "C01"
        # # self.produce.plot.valsGroupBy[2] = "C02"
        # self.produce.plot.valsGroupBy[3] = "C03"

        for key, value in self.produce.plot.valsGroupBy.items():
            self.produce.plot.plotElements[key] = PlotElementConfig(
                plotElement=ScatterElement(),
                title=value,
            )
