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
    "SumAction",
    "MedianHistAction",
    "IqrHistAction",
    "DivideScalar",
    "RmsAction",
)

import logging
import operator
from math import nan
from typing import cast

import numpy as np
from lsst.pex.config import ChoiceField, Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataSchema, Scalar, ScalarAction, Vector
from ...math import nanMax, nanMean, nanMedian, nanMin, nanSigmaMad, nanStd

log = logging.getLogger(__name__)


class ScalarFromVectorAction(ScalarAction):
    """Calculates a statistic from a single vector."""

    vectorKey = Field[str]("Key of Vector to compute statistic from.")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)


class MedianAction(ScalarFromVectorAction):
    """Calculates the median of the given data."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        values = data[self.vectorKey.format(**kwargs)][mask]
        med = nanMedian(values) if len(values) else np.nan

        return med


class MeanAction(ScalarFromVectorAction):
    """Calculates the mean of the given data."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        values = data[self.vectorKey.format(**kwargs)][mask]
        mean = nanMean(values) if len(values) else np.nan

        return mean


class StdevAction(ScalarFromVectorAction):
    """Calculates the standard deviation of the given data."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return nanStd(data[self.vectorKey.format(**kwargs)][mask])


class RmsAction(ScalarFromVectorAction):
    """Calculates the root mean square of the given data (without subtracting
    the mean as in StdevAction)."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        vector = data[self.vectorKey.format(**kwargs)][mask]
        vector = vector[~np.isnan(vector)]

        return np.sqrt(np.mean(vector**2))


class ValueAction(ScalarFromVectorAction):
    """Extracts the first value from a vector."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        return cast(Scalar, float(data[self.vectorKey.format(**kwargs)][0]))


class SigmaMadAction(ScalarFromVectorAction):
    """Calculates the sigma mad of the given data."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return nanSigmaMad(data[self.vectorKey.format(**kwargs)][mask])


class CountAction(ScalarAction):
    """Performs count actions, with threshold-based filtering.
    The operator is specified as a string, for example, "lt", "le", "ge",
    "gt", "ne", and "eq" for the mathematical operations <, <=, >=, >, !=,
    and == respectively. To count non-NaN values, only pass the column name
    as vector key. To count NaN values, pass threshold = nan (from math.nan).
    Optionally to configure from a YAML file, pass "threshold: !!float nan".
    To compute the number of elements with values less than a given threshold,
    use op="le".
    """

    vectorKey = Field[str]("Key of Vector to count")
    op = ChoiceField[str](
        doc="Operator name string.",
        allowed={
            "lt": "less than threshold",
            "le": "less than or equal to threshold",
            "ge": "greater than or equal to threshold",
            "ne": "not equal to a given value",
            "eq": "equal to a given value",
            "gt": "greater than threshold",
        },
        default="ne",
    )
    threshold = Field[float](doc="Threshold to apply.", default=nan)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        arr = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]

        # Count NaNs and non-NaNs
        if self.threshold == nan:
            if self.op == "eq":
                # Count number of NaNs
                result = np.isnan(arr).sum()
                return cast(Scalar, int(result))
            elif self.op == "ne":
                # Count number of non-NaNs
                result = len(arr) - np.isnan(arr).sum()
                return cast(Scalar, int(result))
            else:
                raise ValueError("Invalid operator for counting NaNs.")
        # Count for given threshold ignoring all NaNs
        else:
            result = arr[~np.isnan(arr)]
            result = cast(
                Scalar,
                int(np.sum(getattr(operator, self.op)(result, self.threshold))),
            )
            return result


class CountUniqueAction(ScalarFromVectorAction):
    """Counts the number of unique rows in a given column."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        values = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]
        count = len(np.unique(values))
        return cast(Scalar, count)


class ApproxFloor(ScalarFromVectorAction):
    """Returns the median of the lowest ten values of the sorted input."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(data[self.vectorKey.format(**kwargs)][mask])  # type: ignore
        x = len(value) // 10
        return nanMedian(value[-x:])


class FracThreshold(ScalarFromVectorAction):
    """Compute the fraction of a distribution above or below a threshold.

    The operator is specified as a string, for example,
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
    percent = Field[bool](doc="Express result as percentage", default=False)
    relative_to_median = Field[bool](doc="Calculate threshold relative to the median?", default=False)
    use_absolute_value = Field[bool](
        doc=(
            "Calculate threshold after taking absolute value. If relative_to_median"
            " is true the absolute value will be applied after the median is subtracted"
        ),
        default=False,
    )

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        values = data[self.vectorKey.format(**kwargs)]
        values = values[mask]  # type: ignore
        values = values[np.logical_not(np.isnan(values))]
        n_values = len(values)
        if n_values == 0:
            return np.nan
        threshold = self.threshold
        # If relative_to_median is set, shift the threshold to be median+thresh
        if self.relative_to_median and len(values) > 0:
            offset = nanMedian(values)
            if np.isfinite(offset):
                values -= offset
        if self.use_absolute_value:
            values = np.abs(values)
        result = cast(
            Scalar,
            float(np.sum(getattr(operator, self.op)(values, threshold)) / n_values),  # type: ignore
        )
        if self.percent:
            return 100.0 * result
        else:
            return result


class MaxAction(ScalarFromVectorAction):
    """Returns the maximum of the given data."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return nanMax(data[self.vectorKey.format(**kwargs)][mask])


class MinAction(ScalarFromVectorAction):
    """Returns the minimum of the given data."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        return nanMin(data[self.vectorKey.format(**kwargs)][mask])


class FracInRange(ScalarFromVectorAction):
    """Compute the fraction of a distribution that is between specified
    minimum and maximum values, and is not NaN.
    """

    maximum = Field[float](doc="The maximum value", default=np.nextafter(np.Inf, 0.0))
    minimum = Field[float](doc="The minimum value", default=np.nextafter(-np.Inf, 0.0))
    percent = Field[bool](doc="Express result as percentage", default=False)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
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


class FracNan(ScalarFromVectorAction):
    """Compute the fraction of vector entries that are NaN."""

    percent = Field[bool](doc="Express result as percentage", default=False)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
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


class SumAction(ScalarFromVectorAction):
    """Returns the sum of all values in the column."""

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        arr = cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]
        return cast(Scalar, np.nansum(arr))


class MedianHistAction(ScalarAction):
    """Calculates the median of the given histogram data."""

    histKey = Field[str]("Key of frequency Vector")
    midKey = Field[str]("Key of bin midpoints Vector")

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.histKey, Vector),
            (self.midKey, Vector),
        )

    def histMedian(self, hist, bin_mid):
        """Calculates the median of a histogram with binned values

        Parameters
        ----------
        hist : `numpy.ndarray`
            Frequency array
        bin_mid : `numpy.ndarray`
            Bin midpoints array

        Returns
        -------
        median : `float`
            Median of histogram with binned values
        """
        cumulative_sum = np.cumsum(hist)
        median_index = np.searchsorted(cumulative_sum, cumulative_sum[-1] / 2)
        median = bin_mid[median_index]
        return median

    def __call__(self, data: KeyedData, **kwargs):
        if len(data[self.histKey.format(**kwargs)]) != 0:
            hist = cast(Vector, data[self.histKey.format(**kwargs)])
            bin_mid = cast(Vector, data[self.midKey.format(**kwargs)])
            med = cast(Scalar, float(self.histMedian(hist, bin_mid)))
        else:
            med = np.nan
        return med


class IqrHistAction(ScalarAction):
    """Calculates the interquartile range of the given histogram data."""

    histKey = Field[str]("Key of frequency Vector")
    midKey = Field[str]("Key of bin midpoints Vector")

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.histKey, Vector),
            (self.midKey, Vector),
        )

    def histIqr(self, hist, bin_mid):
        """Calculates the interquartile range of a histogram with binned values

        Parameters
        ----------
        hist : `numpy.ndarray`
            Frequency array
        bin_mid : `numpy.ndarray`
            Bin midpoints array

        Returns
        -------
        iqr : `float`
            Inter-quartile range of histogram with binned values
        """
        cumulative_sum = np.cumsum(hist)
        liqr_index = np.searchsorted(cumulative_sum, cumulative_sum[-1] / 4)
        uiqr_index = np.searchsorted(cumulative_sum, (3 / 4) * cumulative_sum[-1])
        liqr = bin_mid[liqr_index]
        uiqr = bin_mid[uiqr_index]
        iqr = uiqr - liqr
        return iqr

    def __call__(self, data: KeyedData, **kwargs):
        if len(data[self.histKey.format(**kwargs)]) != 0:
            hist = cast(Vector, data[self.histKey.format(**kwargs)])
            bin_mid = cast(Vector, data[self.midKey.format(**kwargs)])
            iqr = cast(Scalar, float(self.histIqr(hist, bin_mid)))
        else:
            iqr = np.nan
        return iqr


class DivideScalar(ScalarAction):
    """Calculate (A/B) for scalars."""

    actionA = ConfigurableActionField[ScalarAction](doc="Action which supplies scalar A")
    actionB = ConfigurableActionField[ScalarAction](doc="Action which supplies scalar B")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()
        yield from self.actionB.getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        """Return the result of A/B.

        Parameters
        ----------
        data : `KeyedData`

        Returns
        -------
        result : `Scalar`
            The result of dividing A by B.
        """
        scalarA = self.actionA(data, **kwargs)
        scalarB = self.actionB(data, **kwargs)
        if scalarB == 0:
            if scalarA == 0:
                log.warning("Both numerator and denominator are zero! Returning NaN.")
                return np.nan
            else:
                log.warning(
                    "Non-zero scalar divided by zero! Returning %sInf." % ("+" if scalarA > 0 else "-")
                )
                return np.sign(scalarA) * np.inf
        else:
            return scalarA / scalarB
