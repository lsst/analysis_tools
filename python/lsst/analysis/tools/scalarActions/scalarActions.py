from __future__ import annotations

from typing import Iterable

import numpy as np
import scipy.stats as sps
from lsst.pex.config import Field

from ..interfaces import NumberType, ScalarAction, Tabular


class MedianAction(ScalarAction):
    key = Field("Key of column to median", dtype=str)

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        return self.key.format(**kwargs)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> NumberType:
        mask = self.getMask(**kwargs)
        return np.nanmedian(table[self.key.format(**kwargs)][mask])  # type: ignore


class SigmaMadAction(ScalarAction):
    key = Field("Key of column to median", dtype=str)

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        return self.key.format(**kwargs)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> NumberType:
        mask = self.getMask(**kwargs)
        return sps.median_abs_deviation(
            table[self.key(**kwargs)][mask],  # type: ignore
            scale="normal",  # type: ignore
            nan_policy="omit",
        )


class CountAction(ScalarAction):
    key = Field("Key of column to median", dtype=str)

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        return self.key.format(**kwargs)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> NumberType:
        mask = self.getMask(**kwargs)
        return len(table[self.key(**kwargs)][mask])  # type: ignore
