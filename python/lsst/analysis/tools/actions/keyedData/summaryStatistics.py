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

__all__ = ("SummaryStatisticAction",)

from lsst.pex.config import Field

from ...interfaces import KeyedData
from ..scalar import CountAction, MedianAction, SigmaMadAction
from .keyedDataActions import KeyedScalars


class SummaryStatisticAction(KeyedScalars):
    vectorKey = Field[str](doc="Column key to compute scalars")

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(vectorKey=self.vectorKey)
        self.scalarActions.sigmaMad = SigmaMadAction(vectorKey=self.vectorKey)
        self.scalarActions.count = CountAction(vectorKey=self.vectorKey)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        return super().__call__(data, **(kwargs | dict(mask=mask)))
