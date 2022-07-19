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

import logging
from typing import cast

import numpy as np
from astropy import units as u
from lsst.pex.config import DictField, Field
from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction
from .selectors import VectorSelector

_LOG = logging.getLogger(__name__)


class DownselectVector(VectorAction):
    """Get a vector from KeyedData, apply specified selector, return the
    shorter Vector.
    """

    vectorKey = Field[str](doc="column key to load from KeyedData")

    selector = ConfigurableActionField(doc="Action which returns a selection mask", default=VectorSelector)

    def getInputSchema(self) -> KeyedDataSchema:
        yield (self.vectorKey, Vector)
        yield from cast(VectorAction, self.selector).getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        mask = cast(VectorAction, self.selector)(data, **kwargs)
        return cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]


class MagColumnNanoJansky(VectorAction):
    vectorKey = Field[str](doc="column key to use for this transformation")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        with np.warnings.catch_warnings():  # type: ignore
            np.warnings.filterwarnings("ignore", r"invalid value encountered")  # type: ignore
            np.warnings.filterwarnings("ignore", r"divide by zero")  # type: ignore
            vec = cast(Vector, data[self.vectorKey.format(**kwargs)])
            return np.array(-2.5 * np.log10((vec * 1e-9) / 3631.0))  # type: ignore


class FractionalDifference(VectorAction):
    """Calculate (A-B)/B"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
    actionB = ConfigurableActionField(doc="Action which supplies vector B", dtype=VectorAction)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema()  # type: ignore
        yield from self.actionB.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return (vecA - vecB) / vecB


class LoadVector(VectorAction):
    """Load and return a Vector from KeyedData"""

    vectorKey = Field[str](doc="Key of vector which should be loaded")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return np.array(cast(Vector, data[self.vectorKey.format(**kwargs)]))


class MagDiff(VectorAction):
    """Calculate the difference between two magnitudes;
    each magnitude is derived from a flux column.
    Parameters
    ----------
    TO DO:
    Returns
    -------
    The magnitude difference in milli mags.
    Notes
    -----
    The flux columns need to be in units (specifiable in
    the fluxUnits1 and 2 config options) that can be converted
    to janskies. This action doesn't have any calibration
    information and assumes that the fluxes are already
    calibrated.
    """

    col1 = Field[str](doc="Column to subtract from")
    fluxUnits1 = Field[str](doc="Units for col1", default="nanojansky")
    col2 = Field[str](doc="Column to subtract")
    fluxUnits2 = Field[str](doc="Units for col2", default="nanojansky")
    returnMillimags = Field[bool](doc="Use millimags or not?", default=True)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.col1, Vector), (self.col2, Vector))

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        flux1 = np.array(data[self.col1.format(**kwargs)]) * u.Unit(self.fluxUnits1)
        mag1 = flux1.to(u.ABmag)

        flux2 = np.array(data[self.col2.format(**kwargs)]) * u.Unit(self.fluxUnits2)
        mag2 = flux2.to(u.ABmag)

        magDiff = mag1 - mag2

        if self.returnMillimags:
            magDiff = magDiff.to(u.mmag)

        return np.array(magDiff.value)


class ExtinctionCorrectedMagDiff(VectorAction):
    """Compute the difference between two magnitudes and correct for extinction
    By default bands are derived from the <band>_ prefix on flux columns,
    per the naming convention in the Object Table:
    e.g. the band of 'g_psfFlux' is 'g'. If column names follow another
    convention, bands can alternatively be supplied via the band1 or band2
    config parameters.
    If band1 and band2 are supplied, the flux column names are ignored.
    """

    magDiff = ConfigurableActionField(
        doc="Action that returns a difference in magnitudes", default=MagDiff, dtype=VectorAction
    )
    ebvCol = Field[str](doc="E(B-V) Column Name", default="ebv")
    band1 = Field[str](
        doc="Optional band for magDiff.col1. Supercedes column name prefix",
        optional=True,
        default=None,
    )
    band2 = Field[str](
        doc="Optional band for magDiff.col2. Supercedes column name prefix",
        optional=True,
        default=None,
    )
    extinctionCoeffs = DictField[str, float](
        doc="Dictionary of extinction coefficients for conversion from E(B-V) to extinction, A_band."
        "Key must be the band",
        optional=True,
        default=None,
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return self.magDiff.getInputSchema() + ((self.ebvCol, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        diff = self.magDiff(data, **kwargs)
        if not self.extinctionCoeffs:
            _LOG.warning("No extinction Coefficients. Not applying extinction correction")
            return diff

        col1Band = self.band1 if self.band1 else self.magDiff.col1.split("_")[0]
        col2Band = self.band2 if self.band2 else self.magDiff.col2.split("_")[0]

        # Return plain MagDiff with warning if either coeff not found
        for band in (col1Band, col2Band):
            if band not in self.extinctionCoeffs:
                _LOG.warning(
                    "%s band not found in coefficients dictionary: %s" " Not applying extinction correction",
                    band,
                    self.extinctionCoeffs,
                )
                return diff

        av1: float = self.extinctionCoeffs[col1Band]
        av2: float = self.extinctionCoeffs[col2Band]

        ebv = data[self.ebvCol]
        # Ignore type until a more complete Vector protocol
        correction = (av1 - av2) * ebv * u.mag  # type: ignore

        if self.magDiff.returnMillimags:
            correction = correction.to(u.mmag)

        return np.array(diff - correction.value)
