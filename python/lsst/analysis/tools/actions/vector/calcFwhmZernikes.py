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

__all__ = ("CalcFwhmZernikesBase", "CalcFwhmZernikesLsst", "CalcFwhmZernikesLatiss")

import numpy as np
from lsst.pex.config import Field
from ...interfaces import KeyedData, KeyedDataSchema, Vector, ScalarAction, Scalar


class CalcFwhmZernikesBase(ScalarAction):

    vectorKey = Field[str](doc="Vector on which to compute statistics")
    # conversion_factors need to be assigned in daughter class
    conversion_factors = None

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:

        results = np.sqrt(np.sum((data[self.vectorKey] * self.conversion_factors) ** 2.0))

        return results


class CalcFwhmZernikesLsst(CalcFwhmZernikesBase):

    conversion_factors = np.array(
        [
            0.751,
            0.271,
            0.271,
            0.819,
            0.819,
            0.396,
            0.396,
            1.679,
            0.937,
            0.937,
            0.517,
            0.517,
            1.755,
            1.755,
            1.089,
            1.089,
            0.635,
            0.635,
            2.810,
        ]
    )


class CalcFwhmZernikesLatiss(CalcFwhmZernikesBase):

    conversion_factors = np.array(
        [
            3.395,
            1.969,
            1.969,
            4.374,
            4.374,
            2.802,
            2.802,
            7.592,
            5.726,
            5.726,
            3.620,
            3.620,
            8.696,
            8.696,
            7.146,
            7.146,
            4.434,
            4.434,
            12.704,
        ]
    )
