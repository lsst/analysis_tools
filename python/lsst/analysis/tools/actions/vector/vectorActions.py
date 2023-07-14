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
    "LoadVector",
    "DownselectVector",
    "MultiCriteriaDownselectVector",
    "ConvertFluxToMag",
    "ConvertUnits",
    "CalcSn",
    "MagDiff",
    "ExtinctionCorrectedMagDiff",
    "PerGroupStatistic",
    "ResidualWithPerGroupStatistic",
    "RAcosDec",
)

import logging
import warnings
from typing import Optional, cast

import numpy as np
import pandas as pd
from astropy import units as u
from lsst.pex.config import DictField, Field
from lsst.pex.config.configurableActions import ConfigurableActionField, ConfigurableActionStructField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction
from .selectors import VectorSelector

_LOG = logging.getLogger(__name__)

# Basic vectorActions


class LoadVector(VectorAction):
    """Load and return a Vector from KeyedData."""

    vectorKey = Field[str](doc="Key of vector which should be loaded")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return np.array(cast(Vector, data[self.vectorKey.format(**kwargs)]))


class DownselectVector(VectorAction):
    """Get a vector from KeyedData, apply specified selector, return the
    shorter Vector.
    """

    vectorKey = Field[str](doc="column key to load from KeyedData")

    selector = ConfigurableActionField[VectorAction](
        doc="Action which returns a selection mask", default=VectorSelector
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield (self.vectorKey, Vector)
        yield from cast(VectorAction, self.selector).getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        mask = cast(VectorAction, self.selector)(data, **kwargs)
        return cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]


class MultiCriteriaDownselectVector(VectorAction):
    """Get a vector from KeyedData, apply specified set of selectors with AND
    logic, and return the shorter Vector.
    """

    vectorKey = Field[str](doc="column key to load from KeyedData")

    selectors = ConfigurableActionStructField[VectorAction](
        doc="Selectors for selecting rows, will be AND together",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield (self.vectorKey, Vector)
        for action in self.selectors:
            yield from action.getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        mask: Optional[Vector] = None
        for selector in self.selectors:
            subMask = selector(data, **kwargs)
            if mask is None:
                mask = subMask
            else:
                mask *= subMask  # type: ignore
        return cast(Vector, data[self.vectorKey.format(**kwargs)])[mask]


# Astronomical vectorActions


class CalcSn(VectorAction):
    """Calculate the signal-to-noise ratio from a single flux vector."""

    fluxType = Field[str](doc="Flux type (vector key) to calculate the S/N.", default="{band}_psfFlux")
    uncertaintySuffix = Field[str](
        doc="Suffix to add to fluxType to specify the uncertainty column", default="Err"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield self.fluxType, Vector
        yield f"{self.fluxType}{self.uncertaintySuffix}", Vector

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        signal = np.array(data[self.fluxType.format(**kwargs)])
        noise = np.array(data[f"{self.fluxType}{self.uncertaintySuffix}".format(**kwargs)])
        sn = signal / noise

        return np.array(sn)


class ConvertFluxToMag(VectorAction):
    """Turn nano janskies into magnitudes."""

    vectorKey = Field[str](doc="Key of flux vector to convert to mags")
    fluxUnit = Field[str](doc="Astropy unit of flux vector", default="nJy")
    returnMillimags = Field[bool](doc="Use millimags or not?", default=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"invalid value encountered")
            warnings.filterwarnings("ignore", r"divide by zero")
            vec = cast(Vector, data[self.vectorKey.format(**kwargs)])
            mags = (np.array(vec) * u.Unit(self.fluxUnit)).to(u.ABmag).value  # type: ignore
            if self.returnMillimags:
                mags *= 1000
            return mags


class ConvertUnits(VectorAction):
    """Convert the units of a vector."""

    buildAction = ConfigurableActionField(doc="Action to build vector", default=LoadVector)
    inUnit = Field[str](doc="input Astropy unit")
    outUnit = Field[str](doc="output Astropy unit")

    def getInputSchema(self) -> KeyedDataSchema:
        return tuple(self.buildAction.getInputSchema())

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        dataWithUnit = self.buildAction(data, **kwargs) * u.Unit(self.inUnit)
        return dataWithUnit.to(self.outUnit).value


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

    magDiff = ConfigurableActionField[VectorAction](
        doc="Action that returns a difference in magnitudes", default=MagDiff
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
            _LOG.debug("No extinction Coefficients. Not applying extinction correction")
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
        correction = np.array((av1 - av2) * ebv) * u.mag  # type: ignore

        if self.magDiff.returnMillimags:
            correction = correction.to(u.mmag)

        return np.array(diff - correction.value)


class RAcosDec(VectorAction):
    """Construct a vector of RA*cos(Dec) in order to have commensurate values
    between RA and Dec."""

    raKey = Field[str](doc="RA coordinate", default="coord_ra")
    decKey = Field[str](doc="Dec coordinate", default="coord_dec")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.decKey, Vector), (self.raKey, Vector))

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        ra = np.array(data[self.raKey])
        dec = np.array(data[self.decKey])
        return ra * np.cos((dec * u.degree).to(u.radian).value)


# Statistical vectorActions


class PerGroupStatistic(VectorAction):
    """Compute per-group statistic values and return result as a vector with
    one element per group. The computed statistic can be any function accepted
    by pandas DataFrameGroupBy.aggregate passed in as a string function name.
    """

    groupKey = Field[str](doc="Column key to use for forming groups", default="obj_index")
    buildAction = ConfigurableActionField[VectorAction](doc="Action to build vector", default=LoadVector)
    func = Field[str](doc="Name of function to be applied per group")

    def getInputSchema(self) -> KeyedDataSchema:
        return tuple(self.buildAction.getInputSchema()) + ((self.groupKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        df = pd.DataFrame({"groupKey": data[self.groupKey], "value": self.buildAction(data, **kwargs)})
        result = df.groupby("groupKey")["value"].aggregate(self.func)
        return np.array(result)


class ResidualWithPerGroupStatistic(VectorAction):
    """Compute residual between individual elements of group and the per-group
    statistic."""

    groupKey = Field[str](doc="Column key to use for forming groups", default="obj_index")
    buildAction = ConfigurableActionField(doc="Action to build vector", default=LoadVector)
    func = Field[str](doc="Name of function to be applied per group", default="mean")

    def getInputSchema(self) -> KeyedDataSchema:
        return tuple(self.buildAction.getInputSchema()) + ((self.groupKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        values = self.buildAction(data, **kwargs)
        df = pd.DataFrame({"groupKey": data[self.groupKey], "value": values})
        result = df.groupby("groupKey")["value"].aggregate(self.func)

        joinedDf = df.join(result, on="groupKey", validate="m:1", lsuffix="_individual", rsuffix="_group")

        result = joinedDf["value_individual"] - joinedDf["value_group"]
        return np.array(result)
