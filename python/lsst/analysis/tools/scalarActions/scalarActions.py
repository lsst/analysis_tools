from __future__ import annotations

import numpy as np
import scipy.stats as sps
from lsst.pex.config import Field

from ..interfaces import Scalar, ScalarAction, KeyedData, KeyedDataSchema, Vector


class MedianAction(ScalarAction):
    colKey = Field("Key of Vector to median", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.colKey.format(**kwargs), Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return np.nanmedian(data[self.colKey.format(**kwargs)][mask])  # type: ignore


class SigmaMadAction(ScalarAction):
    colKey = Field("Key of Vector to median", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.colKey.format(**kwargs), Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return sps.median_abs_deviation(
            data[self.colKey.format(**kwargs)][mask],  # type: ignore
            scale="normal",  # type: ignore
            nan_policy="omit",
        )


class CountAction(ScalarAction):
    colKey = Field("Key of Vector to median", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.colKey.format(**kwargs), Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return len(data[self.colKey.format(**kwargs)][mask])  # type: ignore
