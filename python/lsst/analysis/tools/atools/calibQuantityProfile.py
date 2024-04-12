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

__all__ = (
    "CalibQuantityBaseTool",
    "CalibQuantityAmpProfileScatterTool",
    "CalibQuantityAmpProfileHistTool",
    "CalibAmpScatterTool",
    "CalibDivisaderoScatterTool",
    "CalibPtcCovarScatterTool",
)

from typing import cast

from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ..actions.plot.elements import HistElement, ScatterElement
from ..actions.plot.gridPlot import GridPanelConfig, GridPlot
from ..interfaces import AnalysisTool, KeyedData, KeyedDataAction, PlotElement, Vector


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


class SingleValueRepacker(KeyedDataAction):
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
        uniquePanelKeys = list(set(data[self.panelKey]))

        # Loop over data vector to repack information as it is expected.
        for i in range(len(uniquePanelKeys)):
            repackedData[f"{uniquePanelKeys[i]}_x"] = []
            repackedData[f"{uniquePanelKeys[i]}"] = []

        panelVec = cast(Vector, data[self.panelKey])
        dataVec = cast(Vector, data[self.dataKey])
        quantityVec = cast(Vector, data[self.quantityKey])

        for i in range(len(panelVec)):
            repackedData[f"{panelVec[i]}_x"].append(dataVec[i])
            repackedData[f"{panelVec[i]}"].append(quantityVec[i])

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


class CalibQuantityBaseTool(AnalysisTool):
    parameterizedBand: bool = False

    plotElement = ConfigurableActionField[PlotElement](
        doc="Plot element.",
    )

    def setDefaults(self):
        super().setDefaults()

        # Repack the input data into a usable format
        self.prep = PrepRepacker()

        # A simple pass-through process action to keep the data unchanged
        self.process = PassThrough()

        # Plot the repacked data in a 4x4 grid
        self.produce.plot = GridPlot()
        self.produce.plot.panels = {}
        self.produce.plot.numRows = 4
        self.produce.plot.numCols = 4

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

    def finalize(self):
        super().finalize()

        # Configure each panel
        for key, value in self.produce.plot.valsGroupBy.items():
            gridPanelConfig = GridPanelConfig(
                plotElement=self.plotElement,
                title={"label": str(value), "fontsize": "10"},
                titleY=0.85,
            )
            self.produce.plot.panels[key] = gridPanelConfig


class CalibQuantityAmpProfileScatterTool(CalibQuantityBaseTool):
    def setDefaults(self):
        super().setDefaults()
        self.plotElement = ScatterElement()
        self.prep.panelKey = "amplifier"
        self.prep.dataKey = "mjd"


class CalibQuantityAmpProfileHistTool(CalibQuantityBaseTool):
    def setDefaults(self):
        super().setDefaults()
        self.plotElement = HistElement()
        self.prep.panelKey = "amplifier"
        self.prep.dataKey = "mjd"


class CalibAmpScatterTool(CalibQuantityBaseTool):
    def setDefaults(self):
        super().setDefaults()
        self.plotElement = ScatterElement()

        # Repack the input data into a usable format
        self.prep = SingleValueRepacker()

        self.produce.plot = GridPlot()
        self.produce.plot.panels = {}
        self.produce.plot.numRows = 4
        self.produce.plot.numCols = 4

        # Values to use for x-axis data
        self.produce.plot.xDataKeys = {
            0: "C00_x",
            1: "C01_x",
            2: "C02_x",
            3: "C03_x",
            4: "C04_x",
            5: "C05_x",
            6: "C06_x",
            7: "C07_x",
            8: "C10_x",
            9: "C11_x",
            10: "C12_x",
            11: "C13_x",
            12: "C14_x",
            13: "C15_x",
            14: "C16_x",
            15: "C17_x",
        }

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

        self.prep.panelKey = "amplifier"
        self.prep.dataKey = "mjd"

    def finalize(self):
        super().finalize()

        # Configure each panel
        for key, value in self.produce.plot.valsGroupBy.items():
            gridPanelConfig = GridPanelConfig(
                plotElement=self.plotElement,
                title={"label": str(value), "fontsize": "10"},
                titleY=0.85,
            )
            self.produce.plot.panels[key] = gridPanelConfig
            self.produce.plot.panels[key].plotElement.xKey = f"{value}_x"
            self.produce.plot.panels[key].plotElement.valsKey = value


class CalibDivisaderoScatterTool(CalibQuantityBaseTool):
    def setDefaults(self):
        super().setDefaults()
        self.plotElement = ScatterElement()

        # Repack the input data into a usable format
        self.prep = SingleValueRepacker()

        self.produce.plot = GridPlot()
        self.produce.plot.panels = {}
        self.produce.plot.numRows = 2
        self.produce.plot.numCols = 1

        # Values to group by to distinguish between data in differing panels
        self.produce.plot.valsGroupBy = {0: "bank1", 1: "bank0"}

        self.produce.plot.xDataKeys = {
            0: "bank1_x",
            1: "bank0_x",
        }
        self.prep.panelKey = "amplifier"
        self.prep.dataKey = "mjd"

    def finalize(self):
        super().finalize()

        # Configure each panel
        for key, value in self.produce.plot.valsGroupBy.items():
            gridPanelConfig = GridPanelConfig(
                plotElement=self.plotElement,
                title={"label": str(value), "fontsize": "10"},
                titleY=0.85,
            )
            self.produce.plot.panels[key] = gridPanelConfig
            self.produce.plot.panels[key].plotElement.xKey = f"{value}_x"
            self.produce.plot.panels[key].plotElement.valsKey = value


class CalibPtcCovarScatterTool(CalibQuantityBaseTool):
    def setDefaults(self):
        super().setDefaults()
        self.plotElement = ScatterElement()

        # Repack the input data into a usable format
        self.prep = PrepRepacker()

        self.produce.plot = GridPlot()
        self.produce.plot.panels = {}
        self.produce.plot.numRows = 3
        self.produce.plot.numCols = 2

        # Values to group by to distinguish between data in differing panels
        self.produce.plot.valsGroupBy = {
            0: "PTC_COV_10",
            1: "PTC_COV_20",
            2: "PTC_COV_01",
            3: "PTC_COV_02",
            4: "PTC_COV_11",
            5: None,
        }

        self.produce.plot.xDataKeys = {
            0: "PTC_PTC_RAW_MEANS",
            1: "PTC_PTC_RAW_MEANS",
            2: "PTC_PTC_RAW_MEANS",
            3: "PTC_PTC_RAW_MEANS",
            4: "PTC_PTC_RAW_MEANS",
            5: None,
        }
        self.prep.panelKey = "amplifier"
        self.prep.dataKey = "mjd"

    def finalize(self):
        super().finalize()

        # Configure each panel
        for key, value in self.produce.plot.valsGroupBy.items():
            gridPanelConfig = GridPanelConfig(
                plotElement=self.plotElement,
                title={"label": str(value), "fontsize": "10"},
                titleY=0.85,
            )
            self.produce.plot.panels[key] = gridPanelConfig
            self.produce.plot.panels[key].plotElement.xKey = "PTC_PTC_RAW_MEANS"
            self.produce.plot.panels[key].plotElement.valsKey = f"PTC_{value}"
