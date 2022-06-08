from __future__ import annotations

from typing import cast

import numpy as np
from lsst.pex.config import Field
from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ..interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction


class MagColumnNanoJansky(VectorAction):
    columnKey = Field(doc="column key to use for this transformation", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.columnKey.format(**kwargs), Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        with np.warnings.catch_warnings():  # type: ignore
            np.warnings.filterwarnings("ignore", r"invalid value encountered")  # type: ignore
            np.warnings.filterwarnings("ignore", r"divide by zero")  # type: ignore
            vec = cast(Vector, data[cast(str, self.columnKey).format(**kwargs)])
            return -2.5 * np.log10((vec * 1e-9) / 3631.0)  # type: ignore


class FractionalDifference(VectorAction):
    """Calculate (A-B)/B"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
    actionB = ConfigurableActionField(doc="Action which supplies vector B", dtype=VectorAction)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema(**kwargs)  # type: ignore
        yield from self.actionB.getInputSchema(**kwargs)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return (vecA - vecB) / vecB


class LoadVector(VectorAction):
    """Load and return a Vector from KeyedData"""
    vectorKey = Field(doc="Key of vector which should be loaded", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((cast(str, self.vectorKey).format(**kwargs), Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return cast(Vector, data[cast(str, self.vectorKey).format(**kwargs)])
