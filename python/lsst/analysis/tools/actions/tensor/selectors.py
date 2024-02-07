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
    "PixelMaskSelector",
)

from typing import Optional, Mapping, cast

import numpy as np
from lsst.pex.config.listField import ListField

from ...interfaces import KeyedData, KeyedDataSchema, Tensor, TensorAction

class PixelMaskSelector(TensorAction):
    """Select pixels based on the image pixel mask.
    """

    maskPlaneKeys = ListField[str](doc="Keys of the mask plane dictionary to evaluate", optional=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((key, Tensor) for key in maskPlaneKeys)

    def __call__(self, data: KeyedData, maskPlaneDict: Mapping[str, int] | None = None, **kwargs) -> Tensor:
        """Select on the given mask plane keys
        """

        if maskPlaneDict is None:
            raise ValueError("maskPlaneDict must be supplied")
        if not self.maskPlaneKeys:
            raise RuntimeError("No mask plane keys provided")
        results: Optional[Tensor] = None

        pixelMask = data['mask']
        for key in self.maskPlaneKeys:
            planeBitMask = 2 ** maskPlaneDict[key]
            temp = np.array(pixelMask & planeBitMask) == planeBitMask
            if results is not None:
                results &= temp
            else:
                results = temp

        return cast(Tensor, results)
