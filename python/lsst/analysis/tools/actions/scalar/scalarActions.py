from __future__ import annotations

from typing import cast

import numpy as np
import scipy.stats as sps
from lsst.pex.config import Field

from ...interfaces import KeyedData, KeyedDataSchema, Scalar, ScalarAction, Vector


class MedianAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return np.nanmedian(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])


class MeanAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector from which to calculate mean")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return np.nanmean(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])


class StdevAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector from which to calculate std deviation")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return np.nanstd(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])


class SigmaMadAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return sps.median_abs_deviation(
            cast(Vector, data[self.vectorKey.format(**kwargs)])[mask],
            scale="normal",  # type: ignore
            nan_policy="omit",
        )


class CountAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        arr = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]
        arr = arr[~np.isnan(arr)]
        return len(arr)  # type: ignore


class ApproxFloor(ScalarAction):
    vectorKey = Field[str](doc="Key for the vector to perform action on", optional=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(data[self.vectorKey.format(**kwargs)][mask])  # type: ignore
        x = len(value) // 10
        return np.nanmedian(value[-x:])
