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

__all__ = (#"NumFoundFakesDiaAllMetric",
        #    "FractionFoundFakesDiaAllMetric",
        #    "NumFakesInjectedDiaAllMetric",
           "FractionFoundFakesDiaSnrMetric",
           "FractionFoundFakesDiaMagMetric")
from lsst.pex.config import Field
import numpy as np

from ..actions.scalar import CountAction, FracThreshold
from ..actions.vector import ThresholdSelector, DownselectVector, RangeSelector
from ..interfaces import AnalysisTool


class FractionFoundFakesDiaSnrMetric(AnalysisTool):
    """Calculate the fraction of fake DIA Sources found within the given SNR range
    """

    parameterizedBand: bool = False

    snrMin = Field[float](
        doc="Minimum SNR for fake sources metric calculation.",
        default=0)

    snrMax = Field[float](
        doc="Maximum SNR for fake sources metric calculation.",
        default=np.inf)

    fluxType = Field[str](
        "Flux type for fake sources metric calculation.",
        default="forced_base_PsfFlux_instFlux")

    def setDefaults(self):
        super().setDefaults()

    def finalize(self):
        # There is no need to calculate the SNR as it is already estimated
        # Select the fake sources using their SNR values for the given range and
        # flux type estimation
        self.process.filterActions.fakeSourcesDiaSrcId = DownselectVector(
            vectorKey="diaSourceId",
            selector=RangeSelector(vectorKey=f"{self.fluxType}_SNR", maximum=self.snrMax, minimum=self.snrMin)
        )

        self.process.calculateActions.numTotalFakeSources = CountAction(vectorKey="fakeSourcesDiaSrcId")

        self.process.calculateActions.numFoundFakeSources = CountAction(
            vectorKey="fakeSourcesDiaSrcId", op='gt', threshold=0)

        self.process.calculateActions.fractionFoundFakesDiaAll = FracThreshold(
            vectorKey="fakeSourcesDiaSrcId", op='gt', threshold=0)

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {
            "numTotalFakeSources": "ct",
            "numFoundFakeSources": "ct",
            "fractionFoundFakesDiaAll": "",
        }


class FractionFoundFakesDiaMagMetric(AnalysisTool):
    """Calculate the fraction of fake DIA Sources found within the given magnitude range
    """

    parameterizedBand: bool = False

    magMin = Field[float](
        doc="Minimum magnitude for fake sources metric calculation.",
        default=18)

    magMax = Field[float](
        doc="Maximum magnitude for fake sources metric calculation.",
        default=22)

    def setDefaults(self):
        super().setDefaults()

    def finalize(self):
        # Selecting the fake sources using the truth magnitude values.
        self.process.filterActions.fakeSourcesDiaSrcId = DownselectVector(
            vectorKey="diaSourceId", selector=RangeSelector(
                vectorKey="mag", maximum=self.magMax, minimum=self.magMin)
        )

        self.process.calculateActions.numTotalFakeSources = CountAction(vectorKey="fakeSourcesDiaSrcId")

        self.process.calculateActions.numFoundFakeSources = CountAction(
            vectorKey="fakeSourcesDiaSrcId", op='gt', threshold=0)

        self.process.calculateActions.fractionFoundFakesDiaAll = FracThreshold(
            vectorKey="fakeSourcesDiaSrcId", op='gt', threshold=0)

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {
            "numTotalFakeSources": "ct",
            "numFoundFakeSources": "ct",
            "fractionFoundFakesDiaAll": "",
        }
