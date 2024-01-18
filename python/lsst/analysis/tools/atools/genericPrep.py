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

__all__ = ("CoaddPrep", "VisitPrep", "ExposurePrep")

from ..actions.vector import CoaddPlotFlagSelector, SnSelector, VisitPlotFlagSelector
from ..interfaces import BasePrep, KeyedData, KeyedDataSchema


class CoaddPrep(BasePrep):
    def setDefaults(self):
        super().setDefaults()

        self.selectors.flagSelector = CoaddPlotFlagSelector()
        self.selectors.flagSelector.bands = []

        self.selectors.snSelector = SnSelector()
        self.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.selectors.snSelector.threshold = 100


class VisitPrep(BasePrep):
    def setDefaults(self):
        super().setDefaults()

        self.selectors.flagSelector = VisitPlotFlagSelector()

        self.selectors.snSelector = SnSelector()
        self.selectors.snSelector.fluxType = "psfFlux"
        self.selectors.snSelector.threshold = 100


class ExposurePrep(BasePrep):
    def setDefaults(self):
        super().setDefaults()

    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        pass

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        return data
