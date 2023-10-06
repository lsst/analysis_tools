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
)

import logging
from typing import cast

import numpy as np

from ...interfaces import Vector, VectorAction
from ...statistics import nansigmaMad

_LOG = logging.getLogger(__name__)


class MinMax(VectorAction):

    def __call__(self, data: Vector, **kwargs) -> Vector:
        return cast(Vector, [np.min(data), np.max(data)])


class Med2Mad(VectorAction):

    def __call__(self, data: Vector, **kwargs) -> Vector:
        med = np.nanmedian(data)
        mad = nansigmaMad(data)
        cmin = med - 2 * mad
        cmax = med + 2 * mad
        return cast(Vector, [cmin, cmax])
