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

__all__ = (
    "CalibQuantityBaseTool",
    "CalibQuantityAmpProfileScatterTool",
    "CalibQuantityAmpProfileHistTool",
    "CalibAmpScatterTool",
    "CalibDivisaderoScatterTool",
    "CalibPtcCovarScatterTool",
)

from typing import Optional

import numpy as np
from lsst.pex.config import Field, ListField
from lsst.pex.config.configurableActions import ConfigurableActionField

from ..actions.plot.elements import HistElement, ScatterElement
from ..actions.plot.gridPlot import GridPanelConfig, GridPlot
from ..interfaces import AnalysisTool, KeyedData, KeyedDataAction, KeyedDataSchema, PlotElement, Vector

_CALIB_AMP_NAME_DICT: dict[int, str] = {
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


class PrepRepacker(KeyedDataAction):
    """Prep action to repack data."""

    panelKey = Field[str](
        doc="Panel selector. Data will be separated into multiple panels based on this key.",
    )
    dataKey = Field[str](
        doc="Data selector. Data will be separated into multiple groups in a single panel based on this key.",
    )
    quantityKey = ListField[str](
        doc="Quantity selector. The actual data quantities to be plotted.", minLength=1, optional=False
    )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        repackedData: dict[str, Vector] = {}
        uniquePanelKeys = list(set(data[self.panelKey]))

        for pKey in uniquePanelKeys:
            # Make a boolean array that selects the correct panel data
            sel: np.ndarray = data[self.panelKey] == pKey

            # Setup the x axis
            repackedData[f"{pKey}_x"] = data[self.dataKey][sel]
            for qkey in self.quantityKey:
                # Setup a y axis series for each quantityKey
                repackedData[f"{pKey}_{qkey}"] = data[qkey][sel]
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

    def _get_xKey_dict(self, yKeys: Optional[dict[int, str]] = None) -> dict[int, str]:
        """Generate the dictionary of x axis keys from the y axis ones"""
        if yKeys is None:
            yKeys = self.produce.plot.valsGroupBy
        return {k: f"{v}_x" for k, v in yKeys.items()}

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
        self.produce.plot.valsGroupBy = _CALIB_AMP_NAME_DICT

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
        self.prep = PrepRepacker()

        self.produce.plot = GridPlot()
        self.produce.plot.panels = {}
        self.produce.plot.numRows = 4
        self.produce.plot.numCols = 4

        # Values to use for x-axis data
        # Values to group by to distinguish between data in differing panels
        self.produce.plot.valsGroupBy = _CALIB_AMP_NAME_DICT
        self.produce.plot.xDataKeys = self._get_xKey_dict()

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
        self.prep = PrepRepacker()

        self.produce.plot = GridPlot()
        self.produce.plot.panels = {}
        self.produce.plot.numRows = 2
        self.produce.plot.numCols = 1

        # Values to group by to distinguish between data in differing panels
        self.produce.plot.valsGroupBy = {0: "bank1", 1: "bank0"}

        self.produce.plot.xDataKeys = self._get_xKey_dict()
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
