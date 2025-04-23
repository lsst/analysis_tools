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

__all__ = (
    "cos",
    "divide",
    "isPercent",
    "fluxToMag",
    "nanMax",
    "nanMean",
    "nanMin",
    "nanMedian",
    "nanSigmaMad",
    "nanStd",
    "power",
    "sigmaMad",
    "sin",
    "sqrt",
)

import warnings
from typing import cast

import astropy.units as u
import numpy as np
import scipy.stats as sps

from .interfaces import Scalar, Vector
from .warning_control import (
    filterwarnings_action,
    numpy_all_nan,
    numpy_all_nan_slice,
    numpy_divide_zero_divide,
    numpy_divide_zero_log,
    numpy_divide_zero_log10,
    numpy_dof_zero,
    numpy_invalid_value_cos,
    numpy_invalid_value_divide,
    numpy_invalid_value_log,
    numpy_invalid_value_log10,
    numpy_invalid_value_power,
    numpy_invalid_value_scalar_divide,
    numpy_invalid_value_sin,
    numpy_invalid_value_sqrt,
    numpy_invalid_value_subtract,
    numpy_mean_empty,
)


def cos(values: Scalar | Vector) -> Scalar | Vector:
    """Return the sqrt of values."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_cos)
        result = np.cos(values)
    return result


def divide(dividend: Scalar | Vector, divisor: Scalar | Vector) -> Scalar | Vector:
    """Return dividend/divisor."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_divide_zero_divide)
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_divide)
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_scalar_divide)
        result = dividend / divisor
    return result


def fluxToMag(
    flux: Scalar | Vector,
    flux_unit: u.Unit | str = u.nJy,
    return_millimags: bool = False,
) -> Scalar | Vector:
    """Convert fluxes to magnitudes.

    Parameters
    ----------
    flux
        The flux(es) to convert.
    flux_unit
        The flux unit, as an object or string. Default astropy.units.nJy.
    return_millimags
        Whether to return millimags instead of mags.

    Returns
    -------
    mags
        The magnitude(s) converted from the flux(es).
    """
    if not isinstance(flux_unit, u.Unit):
        flux_unit = u.Unit(flux_unit)
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_divide_zero_log10)
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_log10)
        mag = (np.asarray(flux) * flux_unit).to(u.ABmag).value  # type: ignore
        if return_millimags:
            mag *= 1000
    return mag


def isPercent(value: Scalar) -> bool:
    """Return true if the value is between 0-100"""
    result = 0.0 <= value <= 100.0
    return result


def log(values: Scalar | Vector) -> Scalar | Vector:
    """Return the natural logarithm of values."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_divide_zero_log)
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_log)
        result = np.log(values)
    return result


def log10(values: Scalar | Vector) -> Scalar | Vector:
    """Return the natural logarithm of values."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_divide_zero_log10)
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_log10)
        result = np.log10(values)
    return result


def nanSigmaMad(vector: Vector) -> Scalar:
    """Return the sigma_MAD of a vector."""
    with warnings.catch_warnings():
        # This is needed to catch inf median in sigma_mad
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_subtract)
        result = sps.median_abs_deviation(np.asarray(vector), axis=None, scale="normal", nan_policy="omit")
    return cast(Scalar, result)


def nanMax(vector: Vector) -> Scalar:
    """Return the max of a vector."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_all_nan)
        result = float(np.nanmax(np.asarray(vector)))
    return cast(Scalar, result)


def nanMean(vector: Vector) -> Scalar:
    """Return the mean of a vector."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_all_nan_slice)
        warnings.filterwarnings(filterwarnings_action, numpy_mean_empty)
        result = float(np.nanmean(np.asarray(vector)))
    return cast(Scalar, result)


def nanMedian(vector: Vector) -> Scalar:
    """Return the median of a vector."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_all_nan_slice)
        warnings.filterwarnings(filterwarnings_action, numpy_mean_empty)
        result = float(np.nanmedian(np.asarray(vector)))
    return cast(Scalar, result)


def nanMin(vector: Vector) -> Scalar:
    """Return the max of a vector."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_all_nan)
        result = float(np.nanmin(np.asarray(vector)))
    return cast(Scalar, result)


def nanStd(vector: Vector) -> Scalar:
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_dof_zero)
        result = float(np.nanstd(np.asarray(vector)))
    return cast(Scalar, result)


def power(values: Scalar | Vector, power: float) -> Scalar | Vector:
    """Return the sqrt of values."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_power)
        result = np.power(values, power)
    return result


def sigmaMad(vector: Vector) -> Scalar:
    return cast(Scalar, sps.median_abs_deviation(np.asarray(vector), scale="normal", nan_policy="propagate"))


def sin(values: Scalar | Vector) -> Scalar | Vector:
    """Return the sin of values."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_sin)
        result = np.sin(values)
    return result


def sqrt(values: Scalar | Vector) -> Scalar | Vector:
    """Return the sqrt of values."""
    with warnings.catch_warnings():
        warnings.filterwarnings(filterwarnings_action, numpy_invalid_value_sqrt)
        result = np.sqrt(values)
    return result
