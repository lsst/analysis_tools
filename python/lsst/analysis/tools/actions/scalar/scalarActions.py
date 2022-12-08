from __future__ import annotations

import operator
from typing import cast

import numpy as np
from lsst.pex.config import ChoiceField, Field

from ...interfaces import KeyedData, KeyedDataSchema, Scalar, ScalarAction, Vector
from ...statistics import nansigmaMad


class MedianAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return cast(Scalar, float(np.nanmedian(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))


class MeanAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector from which to calculate mean")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return cast(Scalar, float(np.nanmean(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))


class StdevAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector from which to calculate std deviation")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return cast(Scalar, float(np.nanstd(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))


class SigmaMadAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return cast(
            Scalar,
            float(
                nansigmaMad(
                    data[self.vectorKey.format(**kwargs)][mask],  # type: ignore
                )
            ),
        )


class CountAction(ScalarAction):
    vectorKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        arr = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]
        arr = arr[~np.isnan(arr)]
        return cast(Scalar, len(arr))


class CountUniqueAction(ScalarAction):
    """Counts the number of unique rows in a given column.

    Parameters
    ----------
    data : `KeyedData`

    Returns
    -------
    count : `Scalar`
        The number of unique rows in a given column.
    """

    vectorKey = Field[str](doc="Name of column.")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        values = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]
        count = len(np.unique(values))
        return cast(Scalar, count)


class ApproxFloor(ScalarAction):
    vectorKey = Field[str](doc="Key for the vector to perform action on", optional=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(data[self.vectorKey.format(**kwargs)][mask])  # type: ignore
        x = len(value) // 10
        return cast(Scalar, float(np.nanmedian(value[-x:])))


class FracThreshold(ScalarAction):
    """Compute the fraction of a distribution that is above or below a
    specified threshold. The operator is specified as a string, for example,
    "lt", "le", "ge", "gt" for the mathematical operations <, <=, >=, >. To
    compute the fraction of elements with values less than a given threshold,
    use op="le".
    """

    op = ChoiceField[str](
        doc="Operator name string.",
        allowed={
            "lt": "less than threshold",
            "le": "less than or equal to threshold",
            "ge": "greater than or equal to threshold",
            "gt": "greater than threshold",
        },
    )
    threshold = Field[float](doc="Threshold to apply.")
    vectorKey = Field[str](doc="Name of column")
    percent = Field[bool](doc="Express result as percentage", default=False)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.vectorKey.format(**kwargs), Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        values = data[self.vectorKey.format(**kwargs)]
        values = values[mask]  # type: ignore
        values = values[np.logical_not(np.isnan(values))]
        result = cast(
            Scalar,
            float(np.sum(getattr(operator, self.op)(values, self.threshold)) / len(values)),  # type: ignore
        )
        if self.percent:
            return 100.0 * result
        else:
            return result
