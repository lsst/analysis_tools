from __future__ import annotations

from typing import cast

import numpy as np
from lsst.pex.config import Field

from ..interfaces import Tabular, Vector, VectorAction


class MagColumnNanoJansky(VectorAction):
    columnKey = Field(doc="column key to use for this transformation", dtype=str)

    def __call__(self, table: Tabular, **kwargs) -> Vector:

        with np.warnings.catch_warnings():  # type: ignore
            np.warnings.filterwarnings("ignore", r"invalid value encountered")  # type: ignore
            np.warnings.filterwarnings("ignore", r"divide by zero")  # type: ignore
            return cast(Vector, -2.5 * np.log10((table[self.columnKey] * 1e-9) / 3631.0))
