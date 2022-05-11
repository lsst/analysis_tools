__all__ = ["SNCalculator", "KronFluxDivPsfFlux", "MagDiff", "ColorDiff", "ColorDiffPull",
           "ExtinctionCorrectedMagDiff"]

from lsst.pipe.tasks.configurableActions import ConfigurableActionField
from lsst.pipe.tasks.dataFrameActions import DataFrameAction, DivideColumns, MultiColumnAction
from lsst.pex.config import Field, DictField
from astropy import units as u
import numpy as np
import logging

_LOG = logging.getLogger(__name__)


class SNCalculator(DivideColumns):
    """Calculate the signal to noise by default the i band PSF flux is used"""

    def setDefaults(self):
        super().setDefaults()
        self.colA.column = "i_psfFlux"
        self.colB.column = "i_psfFluxErr"


class KronFluxDivPsfFlux(DivideColumns):
    """Divide the Kron instFlux by the PSF instFlux"""

    def setDefaults(self):
        super().setDefaults()
        self.colA.column = "i_kronFlux"
        self.colB.column = "i_psfFlux"


class MagDiff(MultiColumnAction):
    """Calculate the difference between two magnitudes;
    each magnitude is derived from a flux column.

    Parameters
    ----------
    df : `pandas.core.frame.DataFrame`
        The catalog to calculate the magnitude difference from.

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

    @property
    def columns(self):
        return (self.col1, self.col2)

    def __call__(self, df):
        flux1 = df[self.col1].values * u.Unit(self.fluxUnits1)
        mag1 = flux1.to(u.ABmag)

        flux2 = df[self.col2].values * u.Unit(self.fluxUnits2)
        mag2 = flux2.to(u.ABmag)

        magDiff = mag1 - mag2

        if self.returnMillimags:
            magDiff = magDiff.to(u.mmag)

        return magDiff


class ExtinctionCorrectedMagDiff(DataFrameAction):
    """Compute the difference between two magnitudes and correct for extinction

    By default bands are derived from the <band>_ prefix on flux columns,
    per the naming convention in the Object Table:
    e.g. the band of 'g_psfFlux' is 'g'. If column names follow another
    convention, bands can alternatively be supplied via the band1 or band2
    config parameters.
    If band1 and band2 are supplied, the flux column names are ignored.
    """

    magDiff = ConfigurableActionField(doc="Action that returns a difference in magnitudes",
                                          default=MagDiff, dtype=DataFrameAction)
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


class CalcE(MultiColumnAction):
    """Calculate a complex value representation of the ellipticity

    This is a shape measurement used for doing QA on the ellipticity
    of the sources.

    The complex ellipticity is typically defined as
    E = ((ixx - iyy) + 1j*(2*ixy))/(ixx + iyy) = |E|exp(i*2*theta).

    For plotting purposes we might want to plot |E|*exp(i*theta).
    If `halvePhaseAngle` config parameter is set to `True`, then
    the returned quantity therefore corresponds to |E|*exp(i*theta)
    """

    colXx = Field(doc="The column name to get the xx shape component from.",
                  dtype=str,
                  default="ixx")

    colYy = Field(doc="The column name to get the yy shape component from.",
                  dtype=str,
                  default="iyy")

    colXy = Field(doc="The column name to get the xy shape component from.",
                  dtype=str,
                  default="ixy")

    halvePhaseAngle = Field(doc=("Divide the phase angle by 2? "
                                 "Suitable for quiver plots."),
                            dtype=bool,
                            default=False)

    @property
    def columns(self):
        return (self.colXx, self.colYy, self.colXy)

    def __call__(self, df):
        e = (df[self.colXx] - df[self.colYy]) + 1j*(2*df[self.colXy])
        e /= (df[self.colXx] + df[self.colYy])
        if self.halvePhaseAngle:
            # Ellipiticity is |e|*exp(i*2*theta), but we want to return
            # |e|*exp(i*theta). So we multiply by |e| and take its square root
            # instead of the more expensive trig calls.
            e *= np.abs(e)
            return np.sqrt(e)
        else:
            return e


class CalcEDiff(DataFrameAction):
    """Calculate the difference of two ellipticities as a complex quantity.

    This is a shape measurement used for doing QA on the ellipticity
    of the sources.

    The complex ellipticity difference between E_A and E_B is efined as
    dE = |dE|exp(i*2*theta).

    For plotting purposes we might want to plot |dE|*exp(i*theta).
    If `halvePhaseAngle` config parameter is set to `True`, then
    the returned quantity therefore corresponds to |E|*exp(i*theta)
    """
    colA = ConfigurableActionField(doc="Ellipticity to subtract from",
                                   dtype=MultiColumnAction,
                                   default=CalcE)

    colB = ConfigurableActionField(doc="Ellipticity to subtract",
                                   dtype=MultiColumnAction,
                                   default=CalcE)

    halvePhaseAngle = Field(doc=("Divide the phase angle by 2? "
                                 "Suitable for quiver plots."),
                            dtype=bool,
                            default=False)

    @property
    def columns(self):
        yield from self.colA.columns
        yield from self.colB.columns

    def __call__(self, df):
        eMeas = self.colA(df)
        ePSF = self.colB(df)
        eDiff = eMeas - ePSF
        if self.halvePhaseAngle:
            # Ellipiticity is |e|*exp(i*2*theta), but we want to return
            # |e|*exp(i*theta). So we multiply by |e| and take its square root
            # instead of the more expensive trig calls.
            eDiff *= np.abs(eDiff)
            return np.sqrt(eDiff)
        else:
            return eDiff


class CalcE1(MultiColumnAction):
    """Calculate E1: (ixx - iyy)/(ixx + iyy)
    This is a shape measurement used for doing QA on the ellipticity
    of the sources."""

    colXx = Field(doc="The column name to get the xx shape component from.",
                  dtype=str,
                  default="ixx")

    colYy = Field(doc="The column name to get the yy shape component from.",
                  dtype=str,
                  default="iyy")

    @property
    def columns(self):
        return (self.colXx, self.colYy)

    def __call__(self, df):
        e1 = (df[self.colXx] - df[self.colYy])/(df[self.colXx] + df[self.colYy])

        return e1


class CalcE2(MultiColumnAction):
    """Calculate E2: 2ixy/(ixx+iyy)
    This is a shape measurement used for doing QA on the ellipticity
    of the sources."""

    colXx = Field(doc="The column name to get the xx shape component from.",
                  dtype=str,
                  default="ixx")

    colYy = Field(doc="The column name to get the yy shape component from.",
                  dtype=str,
                  default="iyy")

    colXy = Field(doc="The column name to get the xy shape component from.",
                  dtype=str,
                  default="ixy")

    @property
    def columns(self):
        return (self.colXx, self.colYy, self.colXy)

    def __call__(self, df):
        e2 = 2*df[self.colXy]/(df[self.colXx] + df[self.colYy])
        return e2


class CalcShapeSize(MultiColumnAction):
    """Calculate a size: (ixx*iyy - ixy**2)**0.25

    The square of size measure is typically expressed either as the arithmetic
    mean of the eigenvalues of the moment matrix (trace radius) or as the
    geometric mean of the eigenvalues (determinant radius, computed here).
    Both of these measures give the `sigma^2` parameter for a 2D Gaussian.
    The determinant radius computed here is consistent with the measure
    computed in GalSim:
    http://github.com/GalSim-developers/GalSim/blob/ece3bd32c1ae6ed771f2b489c5ab1b25729e0ea4/galsim/hsm.py#L42

    This is a size measurement used for doing QA on the ellipticity
    of the sources."""

    colXx = Field(doc="The column name to get the xx shape component from.",
                  dtype=str,
                  default="ixx")

    colYy = Field(doc="The column name to get the yy shape component from.",
                  dtype=str,
                  default="iyy")

    colXy = Field(doc="The column name to get the xy shape component from.",
                  dtype=str,
                  default="ixy")

    @property
    def columns(self):
        return (self.colXx, self.colYy, self.colXy)

    def __call__(self, df):
        size = np.power(df[self.colXx]*df[self.colYy] - df[self.colXy]**2, 0.25)
        return size


class ColorDiff(MultiColumnAction):
    """Calculate the difference between two colors;
    each color is derived from two flux columns.

    The color difference is computed as (color1 - color2) with:

    color1 = color1_mag1 - color1_mag2
    color2 = color2_mag1 - color2_mag2

    where color1_mag1 is the magnitude associated with color1_flux1, etc.

    Parameters
    ----------
    df : `pandas.core.frame.DataFrame`
        The catalog to calculate the color difference from.

    Returns
    -------
    The color difference in millimags.

    Notes
    -----
    The flux columns need to be in units that can be converted
    to janskies. This action doesn't have any calibration
    information and assumes that the fluxes are already
    calibrated.
    """
    color1_flux1 = Field(doc="Column for flux1 to determine color1",
                         dtype=str)
    color1_flux1_units = Field(doc="Units for color1_flux1",
                               dtype=str,
                               default="nanojansky")
    color1_flux2 = Field(doc="Column for flux2 to determine color1",
                         dtype=str)
    color1_flux2_units = Field(doc="Units for color1_flux2",
                               dtype=str,
                               default="nanojansky")
    color2_flux1 = Field(doc="Column for flux1 to determine color2",
                         dtype=str)
    color2_flux1_units = Field(doc="Units for color2_flux1",
                               dtype=str,
                               default="nanojansky")
    color2_flux2 = Field(doc="Column for flux2 to determine color2",
                         dtype=str)
    color2_flux2_units = Field(doc="Units for color2_flux2",
                               dtype=str,
                               default="nanojansky")
    return_millimags = Field(doc="Use millimags or not?",
                             dtype=bool,
                             default=True)

    @property
    def columns(self):
        return (self.color1_flux1,
                self.color1_flux2,
                self.color2_flux1,
                self.color2_flux2)

    def __call__(self, df):
        color1_flux1 = df[self.color1_flux1].values*u.Unit(self.color1_flux1_units)
        color1_mag1 = color1_flux1.to(u.ABmag).value

        color1_flux2 = df[self.color1_flux2].values*u.Unit(self.color1_flux2_units)
        color1_mag2 = color1_flux2.to(u.ABmag).value

        color2_flux1 = df[self.color2_flux1].values*u.Unit(self.color2_flux1_units)
        color2_mag1 = color2_flux1.to(u.ABmag).value

        color2_flux2 = df[self.color2_flux2].values*u.Unit(self.color2_flux2_units)
        color2_mag2 = color2_flux2.to(u.ABmag).value

        color1 = color1_mag1 - color1_mag2
        color2 = color2_mag1 - color2_mag2

        color_diff = color1 - color2

        if self.return_millimags:
            color_diff = color_diff*1000

        return color_diff


class ColorDiffPull(ColorDiff):
    """Calculate the difference between two colors, scaled by the color error;
    Each color is derived from two flux columns.

    The color difference is computed as (color1 - color2) with:

    color1 = color1_mag1 - color1_mag2
    color2 = color2_mag1 - color2_mag2

    where color1_mag1 is the magnitude associated with color1_flux1, etc.

    The color difference (color1 - color2) is then scaled by the error on
    the color as computed from color1_flux1_err, color1_flux2_err,
    color2_flux1_err, color2_flux2_err.  The errors on color2 may be omitted
    if the comparison is between an "observed" catalog and a "truth" catalog.

    Parameters
    ----------
    df : `pandas.core.frame.DataFrame`
        The catalog to calculate the color difference from.

    Returns
    -------
    The color difference scaled by the error.

    Notes
    -----
    The flux columns need to be in units that can be converted
    to janskies. This action doesn't have any calibration
    information and assumes that the fluxes are already
    calibrated.
    """
    color1_flux1_err = Field(doc="Error column for flux1 for color1",
                             dtype=str,
                             default="")
    color1_flux2_err = Field(doc="Error column for flux2 for color1",
                             dtype=str,
                             default="")
    color2_flux1_err = Field(doc="Error column for flux1 for color2",
                             dtype=str,
                             default="")
    color2_flux2_err = Field(doc="Error column for flux2 for color2",
                             dtype=str,
                             default="")

    def validate(self):
        super().validate()

        color1_errors = False
        color2_errors = False

        if self.color1_flux1_err and self.color1_flux2_err:
            color1_errors = True
        elif ((self.color1_flux1_err and not self.color1_flux2_err)
              or (not self.color1_flux1_err and self.color1_flux2_err)):
            raise ValueError("Must set both color1_flux1_err and color1_flux2_err if either is set.")
        if self.color2_flux1_err and self.color2_flux2_err:
            color2_errors = True
        elif ((self.color2_flux1_err and not self.color2_flux2_err)
              or (not self.color2_flux1_err and self.color2_flux2_err)):
            raise ValueError("Must set both color2_flux1_err and color2_flux2_err if either is set.")

        if not color1_errors and not color2_errors:
            raise ValueError("Must configure flux errors for at least color1 or color2.")

    @property
    def columns(self):
        columns = (self.color1_flux1,
                   self.color1_flux2,
                   self.color2_flux1,
                   self.color2_flux2)

        if self.color1_flux1_err:
            # Config validation ensures if one is set, both are set.
            columns = columns + (self.color1_flux1_err,
                                 self.color1_flux2_err)

        if self.color2_flux1_err:
            # Config validation ensures if one is set, both are set.
            columns = columns + (self.color2_flux1_err,
                                 self.color2_flux2_err)

        return columns

    def __call__(self, df):
        k = 2.5/np.log(10.)

        color1_flux1 = df[self.color1_flux1].values*u.Unit(self.color1_flux1_units)
        color1_mag1 = color1_flux1.to(u.ABmag).value
        if self.color1_flux1_err:
            color1_mag1_err = k*df[self.color1_flux1_err].values/df[self.color1_flux1].values
        else:
            color1_mag1_err = 0.0

        color1_flux2 = df[self.color1_flux2].values*u.Unit(self.color1_flux2_units)
        color1_mag2 = color1_flux2.to(u.ABmag).value
        if self.color1_flux2_err:
            color1_mag2_err = k*df[self.color1_flux2_err].values/df[self.color1_flux2].values
        else:
            color1_mag2_err = 0.0

        color2_flux1 = df[self.color2_flux1].values*u.Unit(self.color2_flux1_units)
        color2_mag1 = color2_flux1.to(u.ABmag).value
        if self.color2_flux1_err:
            color2_mag1_err = k*df[self.color2_flux1_err].values/df[self.color2_flux1].values
        else:
            color2_mag1_err = 0.0

        color2_flux2 = df[self.color2_flux2].values*u.Unit(self.color2_flux2_units)
        color2_mag2 = color2_flux2.to(u.ABmag).value
        if self.color2_flux2_err:
            color2_mag2_err = k*df[self.color2_flux2_err].values/df[self.color2_flux2].values
        else:
            color2_mag2_err = 0.0

        color1 = color1_mag1 - color1_mag2
        err1_sq = color1_mag1_err**2. + color1_mag2_err**2.
        color2 = color2_mag1 - color2_mag2
        err2_sq = color2_mag1_err**2. + color2_mag2_err**2.

        color_diff = color1 - color2

        pull = color_diff/np.sqrt(err1_sq + err2_sq)

        return pull
