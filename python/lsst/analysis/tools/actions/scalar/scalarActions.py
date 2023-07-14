from __future__ import annotations

__all__ = (
    "MedianAction",
    "MeanAction",
    "StdevAction",
    "ValueAction",
    "SigmaMadAction",
    "CountAction",
    "CountUniqueAction",
    "ApproxFloor",
    "FracThreshold",
    "MaxAction",
    "MinAction",
    "FracInRange",
    "FracNan",
)

import operator
from typing import cast

import numpy as np
from lsst.pex.config import ChoiceField, Field

from ...interfaces import KeyedData, KeyedDataSchema, Scalar, ScalarAction, Vector
from ...statistics import nansigmaMad


class MedianAction(ScalarAction):
    """Calculates the median of the given data."""

    vectorKey = Field[str]("Key of Vector to median")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        if len(data[self.vectorKey.format(**kwargs)][mask]) != 0:
            med = cast(Scalar, float(np.nanmedian(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))
        else:
            med = np.NaN

        return med


class MeanAction(ScalarAction):
    """Calculates the mean of the given data."""

    vectorKey = Field[str]("Key of Vector from which to calculate mean")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        if len(data[self.vectorKey.format(**kwargs)][mask]) != 0:
            mean = cast(Scalar, float(np.nanmean(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))
        else:
            mean = np.NaN

        return mean


class StdevAction(ScalarAction):
    """Calculates the standard deviation of the given data."""

    vectorKey = Field[str]("Key of Vector from which to calculate std deviation")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return cast(Scalar, float(np.nanstd(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))


class ValueAction(ScalarAction):
    """Extracts the first value from a vector."""

    vectorKey = Field[str]("Key of Vector from which to extract the first value")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        return cast(Scalar, float(data[self.vectorKey.format(**kwargs)][0]))


class SigmaMadAction(ScalarAction):
    """Calculates the sigma mad of the given data."""

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
    """Returns the number of non-NaN entries in the given column."""

    vectorKey = Field[str]("Key of Vector to count")

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
    """Returns the median of the lowest ten values of the sorted input."""

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
    relative_to_median = Field[bool](doc="Calculate threshold relative to " "the median?", default=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        values = data[self.vectorKey.format(**kwargs)]
        values = values[mask]  # type: ignore
        values = values[np.logical_not(np.isnan(values))]
        # If relative_to_median is set, shift the threshold to be median+thresh
        if self.relative_to_median:
            threshold = self.threshold + np.median(values)
        else:
            threshold = self.threshold
        result = cast(
            Scalar,
            float(np.sum(getattr(operator, self.op)(values, threshold)) / len(values)),  # type: ignore
        )
        if self.percent:
            return 100.0 * result
        else:
            return result


class MaxAction(ScalarAction):
    """Returns the maximum of the given data."""

    vectorKey = Field[str]("Key of Vector to find maximum")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return cast(Scalar, float(np.max(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))


class MinAction(ScalarAction):
    """Returns the minimum of the given data."""

    vectorKey = Field[str]("Key for the vector to perform action on")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return cast(Scalar, float(np.min(cast(Vector, data[self.vectorKey.format(**kwargs)])[mask])))


class FracInRange(ScalarAction):
    """Compute the fraction of a distribution that is between specified
    minimum and maximum values, and is not NaN.
    """

    vectorKey = Field[str](doc="Name of column")
    maximum = Field[float](doc="The maximum value", default=np.nextafter(np.Inf, 0.0))
    minimum = Field[float](doc="The minimum value", default=np.nextafter(-np.Inf, 0.0))
    percent = Field[bool](doc="Express result as percentage", default=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        """Return the fraction of rows with values within the specified range.

        Parameters
        ----------
        data : `KeyedData`

        Returns
        -------
        result : `Scalar`
            The fraction (or percentage) of rows with values within the
            specified range.
        """
        mask = self.getMask(**kwargs)
        values = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]
        nvalues = len(values)
        values = values[np.logical_not(np.isnan(values))]
        sel_range = (values >= self.minimum) & (values < self.maximum)
        result = cast(
            Scalar,
            float(len(values[sel_range]) / nvalues),  # type: ignore
        )
        if self.percent:
            return 100.0 * result
        else:
            return result


class FracNan(ScalarAction):
    """Compute the fraction of vector entries that are NaN."""

    vectorKey = Field[str](doc="Name of column")
    percent = Field[bool](doc="Express result as percentage", default=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        """Return the fraction of rows with NaN values.

        Parameters
        ----------
        data : `KeyedData`

        Returns
        -------
        result : `Scalar`
            The fraction (or percentage) of rows with NaN values.
        """
        mask = self.getMask(**kwargs)
        values = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]
        nvalues = len(values)
        values = values[np.isnan(values)]
        result = cast(
            Scalar,
            float(len(values) / nvalues),  # type: ignore
        )
        if self.percent:
            return 100.0 * result
        else:
            return result
