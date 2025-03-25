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
    "FractionFoundFakesDiaSnrMetric",
    "FractionFoundFakesDiaMagMetric",
    "FractionFoundFakesAssocDiaSnrMetric",
    "FractionFoundFakesAssocDiaMagMetric",
)

import numpy as np
from lsst.pex.config import Field, ListField

from ..actions.scalar import CountAction, FracThreshold
from ..actions.vector import FlagSelector, MultiCriteriaDownselectVector, RangeSelector
from ..interfaces import AnalysisTool


class FractionFoundFakesDiaSnrMetric(AnalysisTool):
    """Calculate the fraction of fake DIA Sources found within the
    given SNR range"""

    parameterizedBand: bool = False

    snrMin = Field[float](doc="Minimum SNR for fake sources metric calculation.", default=0)

    snrMax = Field[float](doc="Maximum SNR for fake sources metric calculation.", default=np.inf)

    fluxType = Field[str](
        "Flux type for fake sources metric calculation.", default="forced_base_PsfFlux_instFlux"
    )

    fakeFlagsWhenTrue = ListField[str](
        "Flags for fake source cleaning before metrics calculation. Select sources when flags are true.",
        default=[],
    )

    fakeFlagsWhenFalse = ListField[str](
        "Flags for fake source cleaning before metrics calculation. Select sources when flags are false.",
        default=[
            "forced_base_PixelFlags_flag_bad",
            "forced_base_LocalBackground_flag",
            "forced_base_PixelFlags_flag_interpolated",
            "forced_base_PixelFlags_flag_edgeCenter",
        ],
    )

    def setDefaults(self):
        super().setDefaults()

    def finalize(self):
        # There is no need to calculate the SNR as it is already estimated
        # Select the fake sources using their SNR values for the given range
        # and flux type estimation
        self.process.filterActions.fakeSourcesDiaSrcId = MultiCriteriaDownselectVector(
            vectorKey="diaSourceId"
        )

        self.process.filterActions.fakeSourcesDiaSrcId.selectors.snrange = RangeSelector(
            vectorKey=f"{self.fluxType}_SNR", maximum=self.snrMax, minimum=self.snrMin
        )

        self.process.filterActions.fakeSourcesDiaSrcId.selectors.fakeFlags = FlagSelector(
            selectWhenFalse=self.fakeFlagsWhenFalse, selectWhenTrue=self.fakeFlagsWhenTrue
        )

        self.process.calculateActions.numTotalFakeSources = CountAction(vectorKey="fakeSourcesDiaSrcId")

        self.process.calculateActions.numFoundFakeSources = CountAction(
            vectorKey="fakeSourcesDiaSrcId", op="gt", threshold=0
        )

        self.process.calculateActions.fractionFoundFakesDiaAll = FracThreshold(
            vectorKey="fakeSourcesDiaSrcId", op="gt", threshold=0
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {
            "numTotalFakeSources": "ct",
            "numFoundFakeSources": "ct",
            "fractionFoundFakesDiaAll": "",
        }


class FractionFoundFakesDiaMagMetric(AnalysisTool):
    """Calculate the fraction of fake DIA Sources found within the given
    magnitude range"""

    parameterizedBand: bool = False

    magMin = Field[float](doc="Minimum magnitude for fake sources metric calculation.", default=18)

    magMax = Field[float](doc="Maximum magnitude for fake sources metric calculation.", default=22)

    fakeFlagsWhenTrue = ListField[str](
        "Flags for fake source cleaning before metrics calculation.. Select sources when flags are true.",
        default=[],
    )

    fakeFlagsWhenFalse = ListField[str](
        "Flags for fake source cleaning before metrics calculation. Select sources when flags are false.",
        default=[
            "forced_base_PixelFlags_flag_interpolated",
            "forced_base_LocalBackground_flag",
            "forced_base_PixelFlags_flag_bad",
            "forced_base_PixelFlags_flag_edgeCenter",
        ],
    )

    def finalize(self):
        # Selecting the fake sources using the truth magnitude values.
        self.process.filterActions.fakeSourcesDiaSrcId = MultiCriteriaDownselectVector(
            vectorKey="diaSourceId"
        )

        self.process.filterActions.fakeSourcesDiaSrcId.selectors.magrange = RangeSelector(
            vectorKey="mag", maximum=self.magMax, minimum=self.magMin
        )

        self.process.filterActions.fakeSourcesDiaSrcId.selectors.fakeFlags = FlagSelector(
            selectWhenFalse=self.fakeFlagsWhenFalse, selectWhenTrue=self.fakeFlagsWhenTrue
        )

        self.process.calculateActions.numTotalFakeSources = CountAction(vectorKey="fakeSourcesDiaSrcId")

        self.process.calculateActions.numFoundFakeSources = CountAction(
            vectorKey="fakeSourcesDiaSrcId", op="gt", threshold=0
        )

        self.process.calculateActions.fractionFoundFakesDiaAll = FracThreshold(
            vectorKey="fakeSourcesDiaSrcId", op="gt", threshold=0
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {
            "numTotalFakeSources": "ct",
            "numFoundFakeSources": "ct",
            "fractionFoundFakesDiaAll": "",
        }


class FractionFoundFakesAssocDiaSnrMetric(AnalysisTool):
    """Calculate the fraction of fake DIA Sources found within the
    given SNR range"""

    parameterizedBand: bool = False

    snrMin = Field[float](doc="Minimum SNR for fake sources metric calculation.", default=0)

    snrMax = Field[float](doc="Maximum SNR for fake sources metric calculation.", default=np.inf)

    fluxType = Field[str](
        "Flux type for fake sources metric calculation.", default="forced_base_PsfFlux_instFlux"
    )

    fakeFlagsWhenTrue = ListField[str](
        "Flags for fake source cleaning before metrics calculation. Select sources when flags are true.",
        default=[],
    )

    fakeFlagsWhenFalse = ListField[str](
        "Flags for fake source cleaning before metrics calculation. Select sources when flags are false.",
        default=[
            "forced_base_PixelFlags_flag_bad",
            "forced_base_LocalBackground_flag",
            "forced_base_PixelFlags_flag_interpolated",
            "forced_base_PixelFlags_flag_edgeCenter",
        ],
    )

    def setDefaults(self):
        super().setDefaults()

    def finalize(self):
        # There is no need to calculate the SNR as it is already estimated
        # Select the fake sources using their SNR values for the given range
        # and flux type estimation

        self.process.filterActions.fakeSourcesAssocDiaSrcId = MultiCriteriaDownselectVector(
            vectorKey="isAssocDiaSource"
        )

        self.process.filterActions.fakeSourcesAssocDiaSrcId.selectors.snrange = RangeSelector(
            vectorKey=f"{self.fluxType}_SNR", maximum=self.snrMax, minimum=self.snrMin
        )

        self.process.filterActions.fakeSourcesAssocDiaSrcId.selectors.fakeFlags = FlagSelector(
            selectWhenFalse=self.fakeFlagsWhenFalse, selectWhenTrue=self.fakeFlagsWhenTrue
        )

        self.process.calculateActions.numTotalFakeAssocSources = CountAction(
            vectorKey="fakeSourcesAssocDiaSrcId"
        )

        self.process.calculateActions.numFoundFakeAssocSources = CountAction(
            vectorKey="fakeSourcesAssocDiaSrcId", op="gt", threshold=0
        )

        self.process.calculateActions.fractionFoundFakesAssocDiaAll = FracThreshold(
            vectorKey="fakeSourcesAssocDiaSrcId", op="gt", threshold=0
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {
            "numTotalFakeAssocSources": "ct",
            "numFoundFakeAssocSources": "ct",
            "fractionFoundFakesAssocDiaAll": "",
        }


class FractionFoundFakesAssocDiaMagMetric(AnalysisTool):
    """Calculate the fraction of fake sources found in AssocDiasrcs within
    the given magnitude range"""

    parameterizedBand: bool = False

    magMin = Field[float](doc="Minimum magnitude for fake sources metric calculation.", default=18)

    magMax = Field[float](doc="Maximum magnitude for fake sources metric calculation.", default=22)

    fakeFlagsWhenTrue = ListField[str](
        "Flags for fake source cleaning before metrics calculation. Select sources when flags are true.",
        default=[],
    )

    fakeFlagsWhenFalse = ListField[str](
        "Flags for fake source cleaning before metrics calculation. Select sources when flags are false.",
        default=[
            "forced_base_PixelFlags_flag_interpolated",
            "forced_base_LocalBackground_flag",
            "forced_base_PixelFlags_flag_bad",
            "forced_base_PixelFlags_flag_edgeCenter",
        ],
    )

    def finalize(self):
        # Selecting the fake sources using the truth magnitude values.
        self.process.filterActions.fakeSourcesAssocDiaSrcId = MultiCriteriaDownselectVector(
            vectorKey="isAssocDiaSource"
        )

        self.process.filterActions.fakeSourcesAssocDiaSrcId.selectors.magrange = RangeSelector(
            vectorKey="mag", maximum=self.magMax, minimum=self.magMin
        )

        self.process.filterActions.fakeSourcesAssocDiaSrcId.selectors.fakeFlags = FlagSelector(
            selectWhenFalse=self.fakeFlagsWhenFalse, selectWhenTrue=self.fakeFlagsWhenTrue
        )

        self.process.calculateActions.numTotalFakeAssocSources = CountAction(
            vectorKey="fakeSourcesAssocDiaSrcId"
        )

        self.process.calculateActions.numFoundFakeAssocSources = CountAction(
            vectorKey="fakeSourcesAssocDiaSrcId", op="gt", threshold=0
        )

        self.process.calculateActions.fractionFoundFakesAssocDiaAll = FracThreshold(
            vectorKey="fakeSourcesAssocDiaSrcId", op="gt", threshold=0
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {
            "numTotalFakeAssocSources": "ct",
            "numFoundFakeAssocSources": "ct",
            "fractionFoundFakesAssocDiaAll": "",
        }
