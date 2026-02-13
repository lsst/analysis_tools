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
    "CalcRhoStatistics",
    "TreecorrConfig",
)

import logging
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import treecorr  # type: ignore[import]
from deprecated.sphinx import deprecated

from lsst.meas.algorithms.treecorrUtils import TreecorrConfig as TreecorrConfigNew
from lsst.pex.config import ChoiceField, ConfigField, Field

from ...interfaces import KeyedData, KeyedDataAction, Vector
from .calcMomentSize import CalcMomentSize
from .ellipticity import CalcE, CalcEDiff
from .mathActions import FractionalDifference

if TYPE_CHECKING:
    from treecorr import GGCorrelation, KKCorrelation

    from ...interfaces import KeyedDataSchema

_LOG = logging.getLogger(__name__)


# TO DO: Remove TreecorrConfig in here for next major release (DM-47072)
@deprecated(
    reason=(
        "TreecorrConfig is no longer a part of analysis_tools (DM-45899). "
        "Please use lsst.meas.algorithms.treecorrUtils.TreecorrConfig instead."
    ),
    version="v28.0",
    category=FutureWarning,
)
class TreecorrConfig(TreecorrConfigNew):
    pass


class CalcRhoStatistics(KeyedDataAction):
    r"""Calculate rho statistics.

    Rho statistics refer to a collection of correlation functions involving
    PSF ellipticity and size residuals. They quantify the contribution from PSF
    leakage due to errors in PSF modeling to the weak lensing shear correlation
    functions.

    .. _rho_definitions:

    The exact definitions of rho statistics as defined in [1]_ are given below.

    .. math::

       \rho_1(\theta) &= \left\langle
           \delta e^*_{PSF}(x)
           \delta e_{PSF}(x+\theta)
        \right\rangle

       \rho_2(\theta) &= \left\langle
            e^*_{PSF}(x)
            \delta e_{PSF}(x+\theta
        \right\rangle

       \rho_3(\theta) &= \left\langle
            (e^*_{PSF}\frac{\delta T_{PSF}}{T_{PSF}}(x))
            (e_{PSF}\frac{\delta T_{PSF}}{T_{PSF}})(x+\theta)
        \right\rangle

       \rho_4(\theta) &= \left\langle
            \delta e^*_{PSF}(x)
            (e_{PSF}\frac{\delta T_{PSF}}{T_{PSF}})(x+\theta)
        \right\rangle

       \rho_5(\theta) &= \left\langle
            e^*_{PSF}(x)
            (e_{PSF}\frac{\delta T_{PSF}}{T_{PSF}})(x+\theta)
        \right\rangle


    In addition to these five, we also compute the auto-correlation function of
    the fractional size residuals and call it as the :math:`\rho'_3( \theta )`,
    as referred to in Melchior et al. (2015) [2]_.

    .. math::

        \rho'_3(\theta) = \left\langle\frac{\delta T_{PSF}}{T_{PSF}}(x)
                                       \frac{\delta T_{PSF}}{T_{PSF}}(x+\theta)
                           \right\rangle


    The definition of ellipticity used in [1]_ correspond to shear-type,
    which is typically smaller by a factor of 4 than using distortion-type.

    References
    ----------

    .. [1] Jarvis, M., Sheldon, E., Zuntz, J., Kacprzak, T., Bridle, S. L.,
           et. al (2016).
           The DES Science Verification weak lensing shear catalogues
           MNRAS, 460, 2245–2281.
           https://doi.org/10.1093/mnras/stw990;
           https://arxiv.org/abs/1507.05603
    .. [2] Melchior, P., et. al (2015)
           Mass and galaxy distributions of four massive galaxy clusters from
           Dark Energy Survey Science Verification data
           MNRAS, 449, no. 3, pp. 2219–2238.
           https://doi:10.1093/mnras/stv398
           https://arxiv.org/abs/1405.4285
    """

    colRa = Field[str](doc="RA column", default="coord_ra")

    colDec = Field[str](doc="Dec column", default="coord_dec")

    colXx = Field[str](doc="The column name to get the xx shape component from.", default="{band}_ixx")

    colYy = Field[str](doc="The column name to get the yy shape component from.", default="{band}_iyy")

    colXy = Field[str](doc="The column name to get the xy shape component from.", default="{band}_ixy")

    colPsfXx = Field[str](
        doc="The column name to get the PSF xx shape component from.", default="{band}_ixxPSF"
    )

    colPsfYy = Field[str](
        doc="The column name to get the PSF yy shape component from.", default="{band}_iyyPSF"
    )

    colPsfXy = Field[str](
        doc="The column name to get the PSF xy shape component from.", default="{band}_ixyPSF"
    )

    ellipticityType = ChoiceField[str](
        doc="The type of ellipticity to calculate",
        optional=False,
        allowed={
            "distortion": r"Distortion, measured as :math:`(I_{xx}-I_{yy})/(I_{xx}+I_{yy})`",
            "shear": (
                r"Shear, measured as :math:`(I_{xx}-I_{yy})/(I_{xx}+I_{yy}+2\sqrt{I_{xx}I_{yy}-I_{xy}^2})`"
            ),
        },
        default="distortion",
    )

    sizeType = ChoiceField[str](
        doc="The type of size to calculate",
        default="trace",
        allowed={
            "trace": "trace radius",
            "determinant": "determinant radius",
        },
    )

    treecorr = ConfigField[TreecorrConfigNew](
        doc="TreeCorr configuration",
    )

    def setDefaults(self):
        super().setDefaults()
        self.treecorr = TreecorrConfigNew()
        self.treecorr.sep_units = "arcmin"
        self.treecorr.max_sep = 100.0

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.colRa, Vector),
            (self.colDec, Vector),
            (self.colXx, Vector),
            (self.colYy, Vector),
            (self.colXy, Vector),
            (self.colPsfXx, Vector),
            (self.colPsfYy, Vector),
            (self.colPsfXy, Vector),
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        calcEMeas = CalcE(
            colXx=self.colXx,
            colYy=self.colYy,
            colXy=self.colXy,
            ellipticityType=self.ellipticityType,
        )
        calcEpsf = CalcE(
            colXx=self.colPsfXx,
            colYy=self.colPsfYy,
            colXy=self.colPsfXy,
            ellipticityType=self.ellipticityType,
        )

        calcEDiff = CalcEDiff(colA=calcEMeas, colB=calcEpsf)

        calcSizeResidual = FractionalDifference(
            actionA=CalcMomentSize(
                colXx=self.colPsfXx,
                colYy=self.colPsfYy,
                colXy=self.colPsfXy,
                sizeType=self.sizeType,
            ),
            actionB=CalcMomentSize(
                colXx=self.colXx,
                colYy=self.colYy,
                colXy=self.colXy,
                sizeType=self.sizeType,
            ),
        )

        # distortion-type ellipticity has a shear response of 2, so we need to
        # divide by 2 so that the rho-stats do not depend on the
        # ellipticity-type.
        # Note: For distortion, the responsitivity is 2(1 - e^2_{rms}),
        # where e_rms is the root mean square ellipticity per component.
        # This is expected to be small and we ignore it here.
        # This definition of responsitivity is consistent with the definions
        # used in the rho-statistics calculations for the HSC shear catalog
        # papers (Mandelbaum et al. 2018, Li et al., 2022).
        responsitivity = 2.0 if self.ellipticityType == "distortion" else 1.0

        # Call the actions on the data.
        eMEAS = calcEMeas(data, **kwargs)
        if self.ellipticityType == "distortion":
            _LOG.debug("Correction value of responsitivity would be %f", 2 - np.mean(np.abs(eMEAS) ** 2))
        eMEAS /= responsitivity  # type: ignore
        e1, e2 = np.real(eMEAS), np.imag(eMEAS)
        eRes = calcEDiff(data, **kwargs)
        eRes /= responsitivity  # type: ignore
        e1Res, e2Res = np.real(eRes), np.imag(eRes)
        sizeRes = -1 * calcSizeResidual(data, **kwargs)  # sign flip residual to T_psf - T_model

        # Scale the sizeRes by ellipticities
        e1SizeRes = e1 * sizeRes
        e2SizeRes = e2 * sizeRes

        # Package the arguments to capture auto-/cross-correlations for the
        # Rho statistics.
        args = {
            0: (sizeRes, None),
            1: (e1Res, e2Res, None, None),
            2: (e1, e2, e1Res, e2Res),
            3: (e1SizeRes, e2SizeRes, None, None),
            4: (e1Res, e2Res, e1SizeRes, e2SizeRes),
            5: (e1, e2, e1SizeRes, e2SizeRes),
        }

        ra: Vector = data[self.colRa]  # type: ignore
        dec: Vector = data[self.colDec]  # type: ignore

        treecorr_config_dict = self.treecorr.toDict()

        # Swap rng_seed with an rng instance in treecorr config.
        rng = np.random.RandomState(treecorr_config_dict.pop("rng_seed"))
        treecorr_config_dict["rng"] = rng

        # Pass the appropriate arguments to the correlator and build a dict
        rhoStats: Mapping[str, treecorr.BinnedCorr2] = {}
        for rhoIndex in range(1, 6):
            _LOG.info("Calculating rho-%d", rhoIndex)
            rhoStats[f"rho{rhoIndex}"] = self._corrSpin2(  # type: ignore[index]
                ra,
                dec,
                *(args[rhoIndex]),
                treecorr_config_dict=treecorr_config_dict,
            )

        _LOG.info("Calculating rho3alt")
        rhoStats["rho3alt"] = self._corrSpin0(  # type: ignore[index]
            ra,
            dec,
            *(args[0]),
            treecorr_config_dict=treecorr_config_dict,
        )
        return cast(KeyedData, rhoStats)

    @classmethod
    def _corrSpin0(
        cls,
        ra: Vector,
        dec: Vector,
        k1: Vector,
        k2: Vector | None = None,
        raUnits: str = "degrees",
        decUnits: str = "degrees",
        treecorr_config_dict: Mapping[str, Any] | None = None,
    ) -> KKCorrelation:
        """Function to compute correlations between at most two scalar fields.

        This is used to compute rho3alt statistics, given the appropriate
        spin-0 (scalar) fields, usually fractional size residuals.

        Parameters
        ----------
        ra : `numpy.array`
            The right ascension values of entries in the catalog.
        dec : `numpy.array`
            The declination values of entries in the catalog.
        k1 : `numpy.array`
            The primary scalar field.
        k2 : `numpy.array`, optional
            The secondary scalar field.
            Autocorrelation of the primary field is computed if `None`.
        raUnits : `str`, optional
            Unit of the right ascension values. Valid options are
            "degrees", "arcmin", "arcsec", "hours" or "radians".
        decUnits : `str`, optional
            Unit of the declination values. Valid options are
            "degrees", "arcmin", "arcsec", "hours" or "radians".
        treecorr_config_dict: `dict`, optional
            Config dictionary to be passed to `treecorr`
            (`treecorr.KKCorrelation` or `treecorr.Catalog`).

        Returns
        -------
        xy : `treecorr.KKCorrelation`
            A `treecorr.KKCorrelation` object containing the correlation
            function.
        """
        _LOG.debug(
            "No. of entries: %d. The number of pairs in the resulting KKCorrelation cannot exceed %d",
            len(ra),
            len(ra) * (len(ra) - 1) / 2,
        )
        xy = treecorr.KKCorrelation(config=treecorr_config_dict)
        catA = treecorr.Catalog(
            config=treecorr_config_dict,
            ra=ra,
            dec=dec,
            k=k1,
            ra_units=raUnits,
            dec_units=decUnits,
            logger=_LOG,
        )
        if k2 is None:
            # Calculate the auto-correlation
            xy.process(catA)
        else:
            catB = treecorr.Catalog(
                config=treecorr_config_dict,
                ra=ra,
                dec=dec,
                k=k2,
                ra_units=raUnits,
                dec_units=decUnits,
                logger=_LOG,
                patch_centers=catA.patch_centers,
            )
            # Calculate the cross-correlation
            xy.process(catA, catB)

        _LOG.debug("Correlated %d pairs based on the config set.", sum(xy.npairs))
        return xy

    @classmethod
    def _corrSpin2(
        cls,
        ra: Vector,
        dec: Vector,
        g1a: Vector,
        g2a: Vector,
        g1b: Vector | None = None,
        g2b: Vector | None = None,
        raUnits: str = "degrees",
        decUnits: str = "degrees",
        treecorr_config_dict: Mapping[str, Any] | None = None,
    ) -> GGCorrelation:
        """Function to compute correlations between shear-like fields.

        This is used to compute Rho statistics, given the appropriate spin-2
        (shear-like) fields.

        Parameters
        ----------
        ra : `numpy.array`
            The right ascension values of entries in the catalog.
        dec : `numpy.array`
            The declination values of entries in the catalog.
        g1a : `numpy.array`
            The first component of the primary shear-like field.
        g2a : `numpy.array`
            The second component of the primary shear-like field.
        g1b : `numpy.array`, optional
            The first component of the secondary shear-like field.
            Autocorrelation of the primary field is computed if `None`.
        g2b : `numpy.array`, optional
            The second component of the secondary shear-like field.
            Autocorrelation of the primary field is computed if `None`.
        raUnits : `str`, optional
            Unit of the right ascension values. Valid options are
            "degrees", "arcmin", "arcsec", "hours" or "radians".
        decUnits : `str`, optional
            Unit of the declination values. Valid options are
            "degrees", "arcmin", "arcsec", "hours" or "radians".
        treecorr_config_dict : `dict`, optional
            Config dictionary to be passed to `treecorr`
            (`treecorr.GGCorrelation` or `treecorr.Catalog`).

        Returns
        -------
        xy : `treecorr.GGCorrelation`
            A `treecorr.GGCorrelation` object containing the correlation
            function.
        """
        _LOG.debug(
            "No. of entries: %d. The number of pairs in the resulting GGCorrelation cannot exceed %d",
            len(ra),
            len(ra) * (len(ra) - 1) / 2,
        )
        xy = treecorr.GGCorrelation(config=treecorr_config_dict)
        catA = treecorr.Catalog(
            config=treecorr_config_dict,
            ra=ra,
            dec=dec,
            g1=g1a,
            g2=g2a,
            ra_units=raUnits,
            dec_units=decUnits,
            logger=_LOG,
        )
        if g1b is None or g2b is None:
            # Calculate the auto-correlation
            xy.process(catA)
        else:
            catB = treecorr.Catalog(
                config=treecorr_config_dict,
                ra=ra,
                dec=dec,
                g1=g1b,
                g2=g2b,
                ra_units=raUnits,
                dec_units=decUnits,
                logger=_LOG,
                patch_centers=catA.patch_centers,
            )
            # Calculate the cross-correlation
            xy.process(catA, catB)

        _LOG.debug("Correlated %d pairs based on the config set.", sum(xy.npairs))
        return xy
