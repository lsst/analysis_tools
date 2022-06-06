from __future__ import annotations

from typing import cast, Iterable

import numpy as np
from lsst.pex.config import Field
from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ..interfaces import Tabular, Vector, VectorAction


class MagColumnNanoJansky(VectorAction):
    columnKey = Field(doc="column key to use for this transformation", dtype=str)

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        return self.columnKey.format(**kwargs)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> Vector:
        with np.warnings.catch_warnings():  # type: ignore
            np.warnings.filterwarnings("ignore", r"invalid value encountered")  # type: ignore
            np.warnings.filterwarnings("ignore", r"divide by zero")  # type: ignore
            return cast(Vector, -2.5 * np.log10((table[self.columnKey] * 1e-9) / 3631.0))


class FractionalDifference(VectorAction):
    """Calculate (A-B)/B"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
    actionB = ConfigurableActionField(doc="Action which supplies vector B", dtype=VectorAction)

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        yield from self.actionA.getInputColumns(**kwargs)  # type: ignore
        yield from self.actionB.getInputColumns(**kwargs)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> Vector:
        vecA = self.actionA(table, **kwargs)
        vecB = self.actionB(table, **kwargs)
        return (vecA - vecB) / vecB
