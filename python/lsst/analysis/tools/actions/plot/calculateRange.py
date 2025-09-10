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
    "MinMax",
    "Med2Mad",
    "Asinh",
    "Perc",
    "Linear",
)

from typing import cast

import numpy as np
from astropy.visualization import (
    AsinhStretch,
    LinearStretch,
    PercentileInterval,
)
from lsst.pex.config import Field

from ...interfaces import Tensor, TensorAction, Vector, VectorAction
from ...math import nanMax, nanMedian, nanMin, nanSigmaMad


class MinMax(VectorAction):
    """Return the maximum and minimum values of an input vector to use as the
    minimum and maximum values of a colorbar range.

    Parameters
    ----------
    data : `Vector`
        A vector containing the data whose minimum and maximum are to be
        returned.

    Returns
    -------
    A two-element vector containing the minimum and maximum values of `data`.
    """

    def __call__(self, data: Vector, **kwargs) -> Vector:
        if len(data) == 0:
            return cast(Vector, [np.nan, np.nan])
        else:
            return cast(Vector, [nanMin(data), nanMax(data)])


class Med2Mad(VectorAction):
    """Return the median +/- 2*nansigmamad values of an input vector to use
    as the minimum and maximum values of a colorbar range.

    Parameters
    ----------
    data : `Vector`
        A vector containing the data whose median +/- 2*nansigmamad are to
        be returned.

    Returns
    -------
    A two-element vector containing the median +/- 2*nansigmamad values of
    `data`.
    """

    def __call__(self, data: Vector, **kwargs) -> Vector:
        med = nanMedian(data)
        mad = nanSigmaMad(data)
        cmin = med - 2 * mad
        cmax = med + 2 * mad
        return cast(Vector, [cmin, cmax])


class Perc(VectorAction):
    """Return the minimum and maximum values of an input vector after
    excluding a fraction of values from either end of the distribution.

    Parameters
    ----------
    data : `Vector`
        A vector containing the data whose minimum and maximum are to be
        calculated following the exclusion of a fraction of extreme upper
        and lower values.

    Returns
    -------
    A two-element vector containing the minimum and maximum values of
    the input vector, following the exclusion of a fraction of extreme
    upper and lower values.
    """

    percentile = Field[float](
        doc="The fraction of values to keep. The same fraction of values is "
        "eliminated from both ends of the distribution. Default: 97.",
        default=97.0,
    )

    def __call__(self, data, **kwargs):
        return PercentileInterval(self.percentile).get_limits(data)


class Asinh(VectorAction, TensorAction):
    """Transform the input vector/tensor using the asinh stretch.

    Parameters
    ----------
    data : `Vector` | `Tensor`
        A vector or a tensor containing the data to be transformed
        using the asinh stretch.

    Returns
    -------
    A vector or tensor of the same size as the input, transformed
    using the asinh stretch.
    """

    def __call__(self, data: Vector | Tensor, **kwargs) -> Vector | Tensor:
        return AsinhStretch()(data)


class Linear(VectorAction, TensorAction):
    """Transform the input vector/tensor using the linear stretch.

    Parameters
    ----------
    data : `Vector`
        A vector or a tensor containing the data to be transformed
        using the linear stretch.

    Returns
    -------
    A vector or tensor of the same size as the input, transformed
    using the linear stretch.
    """

    intercept = Field[float](doc="The offset of the linear stretch. Default: 0.", default=0.0)

    slope = Field[float](doc="The slope of the linear stretch. Default: 1.", default=1.0)

    def __call__(self, data: Vector | Tensor, **kwargs) -> Vector | Tensor:
        return LinearStretch(self.slope, self.intercept)(data)
