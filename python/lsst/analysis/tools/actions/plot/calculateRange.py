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

__all__ = ("MinMax", "Med2Mad")

from typing import cast

import numpy as np

from ...interfaces import Vector, VectorAction
from ...statistics import nansigmaMad


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
        return cast(Vector, [np.min(data), np.max(data)])


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
        med = np.nanmedian(data)
        mad = nansigmaMad(data)
        cmin = med - 2 * mad
        cmax = med + 2 * mad
        return cast(Vector, [cmin, cmax])
