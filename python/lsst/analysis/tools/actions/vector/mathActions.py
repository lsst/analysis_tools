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
    "ConstantValue",
    "SubtractVector",
    "DivideVector",
    "FractionalDifference",
)

import logging

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction

_LOG = logging.getLogger(__name__)


class ConstantValue(VectorAction):
    """Return a constant scalar value."""

    value = Field[float](doc="A single constant value", optional=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return np.array([self.value])


class SubtractVector(VectorAction):
    """Calculate (A-B)."""

    actionA = ConfigurableActionField[VectorAction](doc="Action which supplies vector A")
    actionB = ConfigurableActionField[VectorAction](doc="Action which supplies vector B")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore
        yield from self.actionB.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return vecA - vecB


class DivideVector(VectorAction):
    """Calculate (A/B)"""

    actionA = ConfigurableActionField[VectorAction](doc="Action which supplies vector A")
    actionB = ConfigurableActionField[VectorAction](doc="Action which supplies vector B")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore
        yield from self.actionB.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return vecA / vecB


class FractionalDifference(VectorAction):
    """Calculate (A-B)/B."""

    actionA = ConfigurableActionField[VectorAction](doc="Action which supplies vector A")
    actionB = ConfigurableActionField[VectorAction](doc="Action which supplies vector B")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore
        yield from self.actionB.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return (vecA - vecB) / vecB
