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

__all__ = ("CPVerifyQuantityProfileTool",)

from typing import cast

from lsst.pex.config import Field

from ..actions.plot.elements.axPlotElement import AxPlotElement
from ..actions.plot.gridPlot import GridPlot, PlotElementConfig
from ..interfaces import AnalysisTool, KeyedData, KeyedDataAction, Vector


class PrepRepacker(KeyedDataAction):
    """Prep action to repack data."""

    panelKey = Field[str](
        doc="Panel selector. Data will be separated into multiple panels based on this key.",
    )
    dataKey = Field[str](
        doc="Data selector. Data will be separated into multiple groups in a single panel based on this key.",
    )
    quantityKey = Field[str](
        doc="Quantity selector. The actual data quantities to be plotted.",
    )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        repackedData = {}
        # Loop over the length of the data vector and repack row by row
        for i in range(len(cast(Vector, data[self.panelKey]))):
            panelVec = cast(Vector, data[self.panelKey])
            dataVec = cast(Vector, data[self.dataKey])
            quantityVec = cast(Vector, data[self.quantityKey])
            repackedData[f"{panelVec[i]}_{dataVec[i]}_{self.quantityKey}"] = quantityVec[i]
        return repackedData

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.panelKey, Vector),
            (self.dataKey, Vector),
            (self.quantityKey, Vector),
        )

    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        pass


class PassThrough(KeyedDataAction):
    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        return data

    def getInputSchema(self) -> KeyedDataSchema:
        # In general this method should be implemented, but here we are ALSO
        # implementing the prep step. We therefore know the method results are
        # not needed because it is doing something special with the inputs.
        return ()


class CPVerifyQuantityProfileTool(AnalysisTool):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        # Repack the input data into a usable format
        self.prep = PrepRepacker()
        self.prep.panelKey = "amplifier"
        self.prep.dataKey = "mjd"
        self.prep.quantityKey = "biasSerialProfile"

        # A simple pass-through process action to keep the data unchanged
        self.process = PassThrough()

        # Plot the repacked data in a 4x4 grid
        self.produce.plot = GridPlot()
        self.produce.plot.plotElements = {}
        self.produce.plot.numRows = 4
        self.produce.plot.numCols = 4
        self.produce.plot.suptitle = {"t": "biasSerialProfile"}

        # Values to group by to distinguish between data in differing panels
        self.produce.plot.valsGroupBy = {
            0: "C00",
            1: "C01",
            2: "C02",
            3: "C03",
            4: "C04",
            5: "C05",
            6: "C06",
            7: "C07",
            8: "C10",
            9: "C11",
            10: "C12",
            11: "C13",
            12: "C14",
            13: "C15",
            14: "C16",
            15: "C17",
        }

        # Set the plot element for each panel to an AxPlotElement
        for key, value in self.produce.plot.valsGroupBy.items():
            plotElementConfig = PlotElementConfig(
                plotElement=AxPlotElement(),
                title={"label": str(value), "fontsize": "10"},
                titleY=0.85,
            )
            self.produce.plot.plotElements[key] = plotElementConfig
