from __future__ import annotations

from typing import cast

from astropy import units as u
import logging
import numpy as np

from lsst.pex.config import Field, DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ..interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction
from .selectors import VectorSelector


_LOG = logging.getLogger(__name__)


class DownselectVector(VectorAction):
    vectorKey = Field(doc="column key to load from KeyedData", dtype=str)

    selector = ConfigurableActionField(
        doc="Action which returns a selection mask",
        default=VectorSelector
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield (cast(str, self.vectorKey).format_map(kwargs), Vector)
        yield from cast(VectorAction, self.selector).getInputSchema(**kwargs)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        mask = cast(VectorAction, self.selector)(data, **kwargs)
        return cast(Vector, data[cast(str, self.vectorKey).format(**kwargs)])[mask]


class MagColumnNanoJansky(VectorAction):
    columnKey = Field(doc="column key to use for this transformation", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.columnKey.format_map(kwargs), Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        with np.warnings.catch_warnings():  # type: ignore
            np.warnings.filterwarnings("ignore", r"invalid value encountered")  # type: ignore
            np.warnings.filterwarnings("ignore", r"divide by zero")  # type: ignore
            vec = cast(Vector, data[cast(str, self.columnKey).format(**kwargs)])
            return -2.5 * np.log10((vec * 1e-9) / 3631.0)  # type: ignore


class FractionalDifference(VectorAction):
    """Calculate (A-B)/B"""

    actionA = ConfigurableActionField(doc="Action which supplies vector A", dtype=VectorAction)
    actionB = ConfigurableActionField(doc="Action which supplies vector B", dtype=VectorAction)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield from self.actionA.getInputSchema(**kwargs)  # type: ignore
        yield from self.actionB.getInputSchema(**kwargs)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        vecA = self.actionA(data, **kwargs)  # type: ignore
        vecB = self.actionB(data, **kwargs)  # type: ignore
        return (vecA - vecB) / vecB


class LoadVector(VectorAction):
    """Load and return a Vector from KeyedData"""
    vectorKey = Field(doc="Key of vector which should be loaded", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((cast(str, self.vectorKey).format_map(kwargs), Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return cast(Vector, data[cast(str, self.vectorKey).format(**kwargs)])


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

    col1 = Field(doc="Column to subtract from", dtype=str)
    fluxUnits1 = Field(doc="Units for col1", dtype=str, default="nanojansky")
    col2 = Field(doc="Column to subtract", dtype=str)
    fluxUnits2 = Field(doc="Units for col2", dtype=str, default="nanojansky")
    returnMillimags = Field(doc="Use millimags or not?", dtype=bool, default=True)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.col1.format_map(kwargs), Vector), (self.col2.format_map(kwargs), Vector))

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        flux1 = data[self.col1.format(**kwargs)] * u.Unit(self.fluxUnits1)
        mag1 = flux1.to(u.ABmag)

        flux2 = data[self.col2.format(**kwargs)] * u.Unit(self.fluxUnits2)
        mag2 = flux2.to(u.ABmag)

        magDiff = mag1 - mag2

        if self.returnMillimags:
            magDiff = magDiff.to(u.mmag)

        return magDiff


class ExtinctionCorrectedMagDiff(VectorAction):
    """Compute the difference between two magnitudes and correct for extinction
    By default bands are derived from the <band>_ prefix on flux columns,
    per the naming convention in the Object Table:
    e.g. the band of 'g_psfFlux' is 'g'. If column names follow another
    convention, bands can alternatively be supplied via the band1 or band2
    config parameters.
    If band1 and band2 are supplied, the flux column names are ignored.
    """

    magDiff = ConfigurableActionField(doc="Action that returns a difference in magnitudes",
                                          default=MagDiff, dtype=VectorAction)
    ebvCol = Field(doc="E(B-V) Column Name", dtype=str, default="ebv")
    band1 = Field(doc="Optional band for magDiff.col1. Supercedes column name prefix",
                  dtype=str, optional=True, default=None)
    band2 = Field(doc="Optional band for magDiff.col2. Supercedes column name prefix",
                  dtype=str, optional=True, default=None)
    extinctionCoeffs = DictField(
        doc="Dictionary of extinction coefficients for conversion from E(B-V) to extinction, A_band."
            "Key must be the band",
        keytype=str, itemtype=float, optional=True,
        default=None)

    @property
    def columns(self):
        return self.magDiff.columns + (self.ebvCol,)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return self.magDiff.getInputSchema(**kwargs) + ((self.ebvCol.format(**kwargs), Vector), )

    def __call__(self, df):
        diff = self.magDiff(df)
        if not self.extinctionCoeffs:
            _LOG.warning("No extinction Coefficients. Not applying extinction correction")
            return diff

        col1Band = self.band1 if self.band1 else self.magDiff.col1.split('_')[0]
        col2Band = self.band2 if self.band2 else self.magDiff.col2.split('_')[0]

        for band in (col1Band, col1Band):
            if band not in self.extinctionCoeffs:
                _LOG.warning("%s band not found in coefficients dictionary: %s"
                             " Not applying extinction correction", band, self.extinctionCoeffs)
                return diff

        av1 = self.extinctionCoeffs[col1Band]
        av2 = self.extinctionCoeffs[col2Band]

        ebv = df[self.ebvCol].values
        correction = (av1 - av2) * ebv * u.mag

        if self.magDiff.returnMillimags:
            correction = correction.to(u.mmag)

        return diff - correction
