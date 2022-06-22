from __future__ import annotations

import numpy as np
import scipy.stats as sps
from lsst.pex.config import Field

from ..interfaces import Scalar, ScalarAction, KeyedData, KeyedDataSchema, Vector


class MedianAction(ScalarAction):
    colKey = Field("Key of Vector to median", dtype=str)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.colKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return np.nanmedian(data[self.colKey.format(**kwargs)][mask])


class SigmaMadAction(ScalarAction):
    colKey = Field("Key of Vector to median", dtype=str)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.colKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return sps.median_abs_deviation(
            data[self.colKey.format(**kwargs)][mask],
            scale="normal",  # type: ignore
            nan_policy="omit",
        )


class CountAction(ScalarAction):
    colKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.colKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return len(data[self.colKey.format(**kwargs)][mask])  # type: ignore


class ApproxFloor(ScalarAction):
    vectorKey = Field[str](doc="Key for the vector to perform action on", optional=False)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.vectorKey.format(**kwargs), Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(data[self.vectorKey.format(**kwargs)][mask])  # type: ignore
        x = int(len(value) / 10)
        return np.nanmedian(value[-x:])
