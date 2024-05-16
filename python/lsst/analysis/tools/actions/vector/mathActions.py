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
    "AddVector",
    "SubtractVector",
    "MultiplyVector",
    "DivideVector",
    "SquareVector",
    "SqrtVector",
    "RaiseFromBaseVector",
    "RaiseToPowerVector",
    "Log10Vector",
    "FractionalDifference",
    "CosVector",
    "SinVector",
)

import logging

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction
from ...math import cos, divide, log10, sin, sqrt

_LOG = logging.getLogger(__name__)


class ConstantValue(VectorAction):
    """Return a constant scalar value."""

    value = Field[float](doc="A single constant value", optional=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return np.array([self.value])


class AddVector(VectorAction):
    """Calculate (A+B)."""

    actionA = ConfigurableActionField[VectorAction](doc="Action which supplies vector A")
    actionB = ConfigurableActionField[VectorAction](doc="Action which supplies vector B")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore
        yield from self.actionB.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return vecA + vecB


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


class MultiplyVector(VectorAction):
    """Calculate (A*B)"""

    actionA = ConfigurableActionField[VectorAction](doc="Action which supplies vector A")
    actionB = ConfigurableActionField[VectorAction](doc="Action which supplies vector B")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore
        yield from self.actionB.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return vecA * vecB


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
        return divide(vecA, vecB)


class SqrtVector(VectorAction):
    """Calculate sqrt(A)"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        return sqrt(vecA)


class SquareVector(VectorAction):
    """Calculate A**2"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        return vecA * vecA


class RaiseFromBaseVector(VectorAction):
    """Calculate n**A"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
    base = Field[float](doc="The base value to raise to the power of vector values")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        return self.base**vecA


class RaiseToPowerVector(VectorAction):
    """Calculate A**n"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
    power = Field[float](doc="The power to raise the vector to")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        return vecA**self.power


class Log10Vector(VectorAction):
    """Calculate log10(A)"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        return log10(vecA)


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
        return divide(vecA - vecB, vecB)


class CosVector(VectorAction):
    """Calculate cos(A)"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        return cos(vecA)


class SinVector(VectorAction):
    """Calculate sin(A)"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        return sin(vecA)
