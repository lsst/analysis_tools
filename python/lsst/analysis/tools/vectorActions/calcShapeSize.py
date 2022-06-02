from __future__ import annotations

from typing import Iterable, cast

import numpy as np
from lsst.pex.config import Field, FieldValidationError
from lsst.pex.config.choiceField import ChoiceField

from ..interfaces import Tabular, Vector, VectorAction


class CalcShapeSize(VectorAction):
    """Calculate a size: (Ixx*Iyy - Ixy**2)**0.25 OR (0.5*(Ixx + Iyy))**0.5
    The square of size measure is typically expressed either as the arithmetic
    mean of the eigenvalues of the moment matrix (trace radius) or as the
    geometric mean of the eigenvalues (determinant radius), which can be
    specified using the ``sizeType`` parameter. Both of these measures give the
    `sigma^2` parameter for a 2D Gaussian.
    Since lensing preserves surface brightness, the determinant radius relates
    the magnification cleanly as it is derived from the area of isophotes, but
    have a slightly higher chance of being NaNs for noisy moment estimates.

    Note
    ----
    This is a size measurement used for doing QA on the ellipticity
    of the sources.
    """

    colXx = Field(
        doc="The column name to get the xx shape component from.",
        dtype=str,
        default="{band}_ixx",
    )

    colYy = Field(
        doc="The column name to get the yy shape component from.",
        dtype=str,
        default="{band}iyy",
    )

    colXy = Field(
        doc="The column name to get the xy shape component from.",
        dtype=str,
        default="{band}_ixy",
        optional=True,
    )

    sizeType = ChoiceField(
        doc="The type of size to calculate",
        dtype=str,
        default="determinant",
        allowed={
            "trace": "trace radius",
            "determinant": "determinant radius",
        },
    )

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        if self.sizeType == "trace":
            return (
                self.colXx.format(**kwargs),  # type: ignore
                self.colYy.format(**kwargs),  # type: ignore
            )
        else:
            return (
                self.colXx.format(**kwargs),  # type: ignore
                self.colYy.format(**kwargs),  # type: ignore
                self.colXy.format(**kwargs),  # type: ignore
            )  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> Vector:
        if self.sizeType == "trace":
            size = np.power(
                0.5
                * (
                    table[self.colXx.format(**kwargs)] + table[self.colYy.format(**kwargs)]  # type: ignore
                ),  # type: ignore
                0.5,
            )
        else:
            size = np.power(
                table[self.colXx.format(**kwargs)] * table[self.colYy.format(**kwargs)]  # type: ignore
                - table[self.colXy] ** 2,
                0.25,
            )

        return cast(Vector, size)

    def validate(self):
        super().validate()
        if self.sizeType == "determinant" and self.colXy is None:
            msg = "colXy is required for determinant-type size"
            raise FieldValidationError(self.__class__.colXy, self, msg)
