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

from typing import cast, TypeVar, Callable, Any, TypeAlias

from lsst.pex.config import Field, ListField
from lsst.pex.config.configurableActions import ConfigurableActionField

from ..actions.plot.elements import HistElement, ScatterElement
from ..actions.plot.gridPlot import GridPanelConfig, GridPlot
from ..interfaces import AnalysisTool, KeyedData, KeyedDataAction, PlotElement, Vector

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

#Dummy class just so we can get the type annotations correct
RepackerLoopFun: TypeAlias = Callable[None, [KeyedData, Any, Any, Any, Any]]

class BaseRepacker(KeyedDataAction):
    """Base class for Data Repacking actions. Essentially Just adds some helper functions"""
    def _repack_loop_helper(self, repackfun: RepackerLoopFun, data: KeyedData) -> KeyedData:
        repackedData: dict[str, Vector] = {}
        quantitiesData = [data[_] for _ in obj.quantityKey]

        for p, d in zip(data[obj.panelKey], data[obj.dataKey]):
            for qName, qD in zip(obj.quantityKey, quantitiesData):
                RepackerLoopFun(repackedData, p, d, qName, qD)
        return repackedData



class PrepRepacker(BaseRepacker):
    """Prep action to repack data."""

    panelKey = Field[str](
        doc="Panel selector. Data will be separated into multiple panels based on this key.",
    )
    dataKey = Field[str](
        doc="Data selector. Data will be separated into multiple groups in a single panel based on this key.",
    )
    quantityKey = ListField[str](
        doc="Quantity selector. The actual data quantities to be plotted.",
        minLength=1, optional=False)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        def rp_loop(rpData: KeyedData, panel, data, quantityName, quantityData):
            qName = f"{panel}_{data}_{quantityName}"
            rpData[qName] = quantityData
        repackedData = self._repack_loop_helper(self, rp_loop, data) 
        return repackedData

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.panelKey, Vector),
            (self.dataKey, Vector),
            (self.quantityKey, Vector),
        )

    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        pass


class SingleValueRepacker(BaseRepacker):
    """Prep action to repack data."""

    panelKey = Field[str](
        doc="Panel selector. Data will be separated into multiple panels based on this key.",
    )
    dataKey = Field[str](
        doc="Data selector. Data will be separated into multiple groups in a single panel based on this key.",
    )
    quantityKey = ListField[str](
        doc="Quantity selector. The actual data quantities to be plotted.",
        minLength=1, optional=False)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        repackedData: dict[str, Vector] = {}

        def rp_loop(rpData: KeyedData, panel, data, quantityName, quantityData):
            if (xlab := f"{panel}_x") not in rpData:
                rpData[xlab] = [data]
            else:
                rpData[xlab].append(data)

            if (lab := f"{panel}_{quantityName}") not in rpData:
                rpData[lab] = [quantityData]
            else:
                rpData[lab].append(quantityData)

        repackedData = self._repack_loop_helper(self, rp_loop, data)
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
        self.prep = SingleValueRepacker()

        self.produce.plot = GridPlot()
        self.produce.plot.panels = {}
        self.produce.plot.numRows = 4
        self.produce.plot.numCols = 4

        # Values to use for x-axis data
        self.produce.plot.xDataKeys = {k : f"{v}_x" for k,v in _CALIB_AMP_NAME_DICT.items()}

        # Values to group by to distinguish between data in differing panels
        self.produce.plot.valsGroupBy = _CALIB_AMP_NAME_DICT
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
