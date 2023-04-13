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
    "CalcE",
    "CalcEDiff",
    "CalcE1",
    "CalcE2",
)

import numpy as np
from lsst.pex.config import ChoiceField, Field, FieldValidationError
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction


class CalcE(VectorAction):
    """Calculate a complex value representation of the ellipticity.

    The complex ellipticity is typically defined as
    :math:`e = |e|\exp{(j*2*theta)} = \\frac{((Ixx - Iyy) + j*(2*Ixy))}{(Ixx + Iyy)}`, where j is
    the square root of -1 and Ixx, Iyy, Ixy are second-order central moments.
    This is sometimes referred to as distortion, and denoted by e = (e1, e2)
    in GalSim (see Eq. 4.4. of Bartelmann and Schneider, 2001).
    The other definition differs in normalization.
    It is referred to as shear, and denoted by g = (g1, g2)
    in GalSim (see Eq. 4.10 of Bartelmann and Schneider (2001). It is defined
    as :math:`g = \\frac{((Ixx - Iyy) + j*(2*Ixy))}{(Ixx + Iyy + 2\sqrt{(Ixx*Iyy - Ixy^{2})})}`.

    The shear measure is unbiased in weak-lensing shear, but may exclude some
    objects in the presence of noisy moment estimates. The distortion measure
    is biased in weak-lensing distortion, but does not suffer from selection
    artifacts.

    References
    ----------
    [1] Bartelmann, M. and Schneider, P., “Weak gravitational lensing”,
    Physics Reports, vol. 340, no. 4–5, pp. 291–472, 2001.
    doi:10.1016/S0370-1573(00)00082-X; https://arxiv.org/abs/astro-ph/9912508

    Notes
    -----

    1. This is a shape measurement used for doing QA on the ellipticity
    of the sources.

    2. For plotting purposes we might want to plot :math:`|E|*\exp{(i*theta)}`.
    If `halvePhaseAngle` config parameter is set to `True`, then
    the returned quantity therefore corresponds to :math:`|E|*\exp{(i*theta)}`.

    See Also
    --------
    CalcE1
    CalcE2
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
    )

    ellipticityType = ChoiceField[str](
        doc="The type of ellipticity to calculate",
        allowed={
            "distortion": ("Distortion, defined as (Ixx - Iyy + 2j*Ixy)/" "(Ixx + Iyy)"),
            "shear": ("Shear, defined as (Ixx - Iyy + 2j*Ixy)/" "(Ixx + Iyy + 2*sqrt(Ixx*Iyy - Ixy**2))"),
        },
        default="distortion",
    )

    halvePhaseAngle = Field[bool](
        doc="Divide the phase angle by 2? Suitable for quiver plots.",
        default=False,
    )

    component = ChoiceField[str](
        doc="Which component of the ellipticity to return. If `None`, return complex ellipticity values.",
        optional=True,
        allowed={
            "1": "e1 or g1 (depending on `ellipticityType`)",
            "2": "e2 or g2 (depending on `ellipticityType`)",
        },
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.colXx, Vector), (self.colXy, Vector), (self.colYy, Vector))

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        e = (data[self.colXx.format(**kwargs)] - data[self.colYy.format(**kwargs)]) + 1j * (
            2 * data[self.colXy.format(**kwargs)]
        )
        denom = data[self.colXx.format(**kwargs)] + data[self.colYy.format(**kwargs)]

        if self.ellipticityType == "shear":
            denom += 2 * np.sqrt(
                data[self.colXx.format(**kwargs)] * data[self.colYy.format(**kwargs)]
                - data[self.colXy.format(**kwargs)] ** 2
            )

        e /= denom

        if self.halvePhaseAngle:
            # Ellipiticity is |e|*exp(i*2*theta), but we want to return
            # |e|*exp(i*theta). So we multiply by |e| and take its square root
            # instead of the more expensive trig calls.
            e *= np.abs(e)
            e = np.sqrt(e)

        if self.component == "1":
            return np.real(e)
        elif self.component == "2":
            return np.imag(e)
        else:
            return e


class CalcEDiff(VectorAction):
    """Calculate the difference of two ellipticities as a complex quantity.

    The complex ellipticity difference between e_A and e_B is defined as
    :math:`e_{A} - e_{B} = de = |de|\exp{(j*2*theta)}`.

    See Also
    --------
    CalcE

    Notes
    -----

    1. This is a shape measurement used for doing QA on the ellipticity
    of the sources.

    2. For plotting purposes we might want to plot :math:`|de|*\exp{(j*theta)}`.
    If `halvePhaseAngle` config parameter is set to `True`, then
    the returned quantity therefore corresponds to :math:`|e|*\exp{(j*theta)}`.
    """

    colA = ConfigurableActionField(
        doc="Ellipticity to subtract from",
        dtype=VectorAction,
        default=CalcE,
    )

    colB = ConfigurableActionField(
        doc="Ellipticity to subtract",
        dtype=VectorAction,
        default=CalcE,
    )

    halvePhaseAngle = Field[bool](
        doc="Divide the phase angle by 2? Suitable for quiver plots.",
        default=False,
    )

    component = ChoiceField[str](
        doc="Which component of the ellipticity to return. If `None`, return complex ellipticity values.",
        optional=True,
        allowed={
            "1": "e1 or g1 (depending on the `ellipiticyType`)",
            "2": "e2 or g2 (depending on the `ellipiticyType`)",
        },
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.colA.getInputSchema()
        yield from self.colB.getInputSchema()

    def validate(self):
        super().validate()
        if self.colA.ellipticityType != self.colB.ellipticityType:
            msg = "Both the ellipticities in CalcEDiff must have the same type."
            raise FieldValidationError(self.colB.__class__.ellipticityType, self, msg)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        eMeas = self.colA(data, **kwargs)
        ePSF = self.colB(data, **kwargs)
        eDiff = eMeas - ePSF
        if self.halvePhaseAngle:
            # Ellipiticity is |e|*exp(i*2*theta), but we want to return
            # |e|*exp(j*theta). So we multiply by |e| and take its square root
            # instead of the more expensive trig calls.
            eDiff *= np.abs(eDiff)
            eDiff = np.sqrt(eDiff)

        if self.component == "1":
            return np.real(eDiff)
        elif self.component == "2":
            return np.imag(eDiff)
        else:
            return eDiff


class CalcE1(VectorAction):
    """Calculate distortion-type :math:`e1 = (Ixx - Iyy)/(Ixx + Iyy)` or
    shear-type :math:`g1 = (Ixx - Iyy)/(Ixx + Iyy + 2sqrt(Ixx*Iyy - Ixy^{2}))`.

    See Also
    --------
    CalcE
    CalcE2

    Note
    ----
    This is a shape measurement used for doing QA on the ellipticity
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

    ellipticityType = ChoiceField[str](
        doc="The type of ellipticity to calculate",
        allowed={
            "distortion": "Distortion, measured as (Ixx - Iyy)/(Ixx + Iyy)",
            "shear": ("Shear, measured as (Ixx - Iyy)/" "(Ixx + Iyy + 2*sqrt(Ixx*Iyy - Ixy**2))"),
        },
        default="distortion",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        if self.ellipticityType == "distortion":
            return (
                (self.colXx, Vector),
                (self.colYy, Vector),
            )
        else:
            return (
                (self.colXx, Vector),
                (self.colYy, Vector),
                (self.colXy, Vector),
            )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        denom = data[self.colXx.format(**kwargs)] + data[self.colYy.format(**kwargs)]
        if self.ellipticityType == "shear":
            denom += 2 * np.sqrt(
                data[self.colXx.format(**kwargs)] * data[self.colYy.format(**kwargs)]
                - data[self.colXy.format(**kwargs)] ** 2
            )
        e1 = (data[self.colXx.format(**kwargs)] - data[self.colYy.format(**kwargs)]) / denom

        return e1

    def validate(self):
        super().validate()
        if self.ellipticityType == "shear" and self.colXy is None:
            msg = "colXy is required for shear-type shear ellipticity"
            raise FieldValidationError(self.__class__.colXy, self, msg)


class CalcE2(VectorAction):
    """Calculate distortion-type :math:`e2 = 2Ixy/(Ixx+Iyy)` or
    shear-type :math:`g2 = 2Ixy/(Ixx+Iyy+2\sqrt(Ixx*Iyy - Ixy^{2}))`.

    See Also
    --------
    CalcE
    CalcE1

    Note
    ----
    This is a shape measurement used for doing QA on the ellipticity
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
    )

    ellipticityType = ChoiceField[str](
        doc="The type of ellipticity to calculate",
        allowed={
            "distortion": "Distortion, defined as 2*Ixy/(Ixx + Iyy)",
            "shear": ("Shear, defined as 2*Ixy/" "(Ixx + Iyy + 2*sqrt(Ixx*Iyy - Ixy**2))"),
        },
        default="distortion",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.colXx, Vector),
            (self.colYy, Vector),
            (self.colXy, Vector),
        )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        denom = data[self.colXx.format(**kwargs)] + data[self.colYy.format(**kwargs)]
        if self.ellipticityType == "shear":
            denom += 2 * np.sqrt(
                data[self.colXx.format(**kwargs)] * data[self.colYy.format(**kwargs)]
                - data[self.colXy.format(**kwargs)] ** 2
            )
        e2 = 2 * data[self.colXy.format(**kwargs)] / denom
        return e2
