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

__all__ = ("CalcShapeSize",)

import numpy as np
from lsst.pex.config import Field, FieldValidationError
from lsst.pex.config.choiceField import ChoiceField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction


class CalcShapeSize(VectorAction):
    r"""Calculate a size: :math:`(I_{xx}I_{yy}-I_{xy}^2)^{\frac{1}{4}}`
    (determinant radius) or :math:`\sqrt{(I_{xx}+I_{yy})/2}`
    (trace radius).

    The square of size measure is typically expressed either as the arithmetic
    mean of the eigenvalues of the moment matrix (trace radius) or as the
    geometric mean of the eigenvalues (determinant radius), which can be
    specified using the `sizeType` parameter. Both of these measures give the
    :math:`\sigma^2` parameter for a 2D Gaussian.

    Since lensing preserves surface brightness, the determinant radius relates
    the magnification cleanly as it is derived from the area of isophotes, but
    have a slightly higher chance of being NaNs for noisy moment estimates.

    Notes
    -----
    This is a size measurement used for doing QA on the ellipticity
    of the sources.
    """

    colXx = Field[str](
        doc="The column name to get the xx shape component from.",
        default="{band}_ixx",
    )

    colYy = Field[str](
        doc="The column name to get the yy shape component from.",
        default="{band}_iyy",
    )

    colXy = Field[str](
        doc="The column name to get the xy shape component from.",
        default="{band}_ixy",
        optional=True,
    )

    sizeType = ChoiceField[str](
        doc="The type of size to calculate",
        default="determinant",
        optional=False,
        allowed={
            "trace": r"Trace radius :math:`\sqrt{(I_{xx}+I_{yy})/2}`",
            "determinant": r"Determinant radius :math:`(I_{xx}I_{yy}-I_{xy}^2)^{\frac{1}{4}}`",
        },
    )

    def getInputSchema(self) -> KeyedDataSchema:
        if self.sizeType == "trace":
            return (
                (self.colXx, Vector),
                (self.colYy, Vector),
            )
        else:
            return (
                (self.colXx, Vector),
                (self.colYy, Vector),
                (self.colXy, Vector),
            )  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        if self.sizeType == "trace":
            size = np.power(
                0.5
                * (
                    data[self.colXx.format(**kwargs)] + data[self.colYy.format(**kwargs)]  # type: ignore
                ),  # type: ignore
                0.5,
            )
        else:
            size = np.power(
                data[self.colXx.format(**kwargs)] * data[self.colYy.format(**kwargs)]  # type: ignore
                - data[self.colXy.format(**kwargs)] ** 2,  # type: ignore
                0.25,
            )

        return size

    def validate(self):
        super().validate()
        if self.sizeType == "determinant" and self.colXy is None:
            msg = "colXy is required for determinant-type size"
            raise FieldValidationError(self.__class__.colXy, self, msg)
