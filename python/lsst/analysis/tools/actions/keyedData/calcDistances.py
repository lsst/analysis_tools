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
__all__ = ("CalcRelativeDistances",)

import astropy.units as u
import esutil
import numpy as np
from smatch import Matcher

from lsst.pex.config import Field

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Vector


class CalcRelativeDistances(KeyedDataAction):
    """Calculate relative distances in a matched catalog.

    Given a catalog of matched sources from multiple visits, this finds all
    pairs of objects at a given separation, then calculates the separation of
    their component source measurements from the individual visits. The RMS of
    these is used to calculate the astrometric relative repeatability metric,
    AMx, while the overall distribution of separations is used to compute the
    ADx and AFx metrics.
    """

    groupKey = Field[str](doc="Column key to use for forming groups", default="obj_index")
    visitKey = Field[str](doc="Column key to use for matching visits", default="visit")
    raKey = Field[str](doc="RA column key", default="coord_ra")
    decKey = Field[str](doc="Dec column key", default="coord_dec")
    annulus = Field[float](doc="Radial distance of the annulus in arcmin", default=5.0)
    width = Field[float](doc="Width of annulus in arcmin", default=2.0)
    threshAD = Field[float](doc="Threshold in mas for AFx calculation.", default=20.0)
    threshAF = Field[float](
        doc="Percentile of differences that can vary by more than threshAD.", default=10.0
    )
    maxPairs = Field[int](
        doc="Maximum number of pairs to use; downsample otherwise.",
        default=100_000,
    )
    bruteForceRadiusThreshold = Field[float](
        doc="Matching radius at which to switch to brute-force matching (deg).",
        default=0.1,
    )
    bruteForceMaxMatch = Field[int](
        doc="Number of objects to sample for brute-force matching.",
        default=10_000,
    )
    randomSeed = Field[int](
        doc="Random seed to use when downsampling.",
        default=12345,
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.groupKey, Vector),
            (self.raKey, Vector),
            (self.decKey, Vector),
            (self.visitKey, Vector),
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        """Run the calculation.

        Parameters
        ----------
        data : KeyedData
            Catalog of data including coordinate, visit, and object group
            information.

        Returns
        -------
        distanceParams : `dict`
            Dictionary of the calculated arrays and metrics with the following
            keys:

            - ``rmsDistances`` : Per-object rms of separations (`np.array`).
            - ``separationResiduals`` : All separations minus per-object median
                (`np.array`)
            - ``AMx`` : AMx metric (`float`).
            - ``ADx`` : ADx metric (`float`).
            - ``AFx`` : AFx metric (`float`).
        """
        distanceParams = {
            "rmsDistances": np.array([]),
            "separationResiduals": np.array([]),
            "AMx": np.nan,
            "ADx": np.nan,
            "AFx": np.nan,
        }

        if len(data[self.groupKey]) == 0:
            return distanceParams

        rng = np.random.RandomState(seed=self.randomSeed)

        def _compressArray(arrayIn):
            h, rev = esutil.stat.histogram(arrayIn, rev=True)
            arrayOut = np.zeros(len(arrayIn), dtype=np.int32)
            (good,) = np.where(h > 0)
            for counter, ind in enumerate(good):
                arrayOut[rev[rev[ind] : rev[ind + 1]]] = counter
            return arrayOut

        groupId = _compressArray(data[self.groupKey])

        nObj = groupId.max() + 1

        # Compute the meanRa/meanDec.
        meanRa = np.zeros(nObj)
        meanDec = np.zeros_like(meanRa)
        nObs = np.zeros_like(meanRa, dtype=np.int64)

        # Check if tract is overlapping ra=0 and rotate if so.
        # We assume a tract is smaller than 60x60 degrees.
        rotation = 0.0
        if np.max(data[self.raKey]) > 330.0 and np.min(data[self.raKey]) < 30.0:
            rotation = 180.0
            raRotated = np.array(data[self.raKey]) - rotation
        else:
            raRotated = np.array(data[self.raKey])

        np.add.at(meanRa, groupId, raRotated)
        np.add.at(meanDec, groupId, np.array(data[self.decKey]))
        np.add.at(nObs, groupId, 1)

        meanRa /= nObs
        meanDec /= nObs
        meanRa += rotation

        D = (self.annulus * u.arcmin).to_value(u.degree)
        width = (self.width * u.arcmin).to_value(u.degree)
        annulus = D + (width / 2) * np.array([-1, +1])

        if annulus[1] >= self.bruteForceRadiusThreshold:
            # Use brute-force matching.
            if len(meanRa) > self.bruteForceMaxMatch:
                bruteForceSelection = rng.choice(len(meanRa), size=self.bruteForceMaxMatch, replace=False)
            else:
                bruteForceSelection = np.arange(len(meanRa))

            i1 = np.tile(bruteForceSelection, len(bruteForceSelection))
            i2 = np.repeat(bruteForceSelection, len(bruteForceSelection))

            d = np.rad2deg(
                sphDist(
                    np.deg2rad(meanRa)[i1],
                    np.deg2rad(meanDec)[i1],
                    np.deg2rad(meanRa)[i2],
                    np.deg2rad(meanDec)[i2],
                ),
            )
        else:
            with Matcher(meanRa, meanDec) as m:
                idx, i1, i2, d = m.query_self(annulus[1], return_indices=True)

        inAnnulus = (d > annulus[0]) & (d < annulus[1])
        i1 = i1[inAnnulus]
        i2 = i2[inAnnulus]

        if len(i1) == 0:
            return distanceParams

        if len(i1) > self.maxPairs:
            # Downsample the pairs.
            selection = rng.choice(len(i1), size=self.maxPairs, replace=False)
            i1 = i1[selection]
            i2 = i2[selection]

        # Match groups and get indices.
        h, rev = esutil.stat.histogram(groupId, rev=True)

        # Match together pairs that have the same visit.
        # It is unfortunate that this requires a loop, but it is not slow.
        # After this matching we have a set of matchedObsInd1/matchedObsInd2
        # that are all individual observations that are in the annulus and
        # share a visit. The matchedPairInd groups all the paired observations
        # of a given pair.
        matchedObsInd1 = []
        matchedObsInd2 = []
        matchedPairInd = []
        for ind in range(len(i1)):
            objInd1 = i1[ind]
            objInd2 = i2[ind]
            obsInd1 = rev[rev[objInd1] : rev[objInd1 + 1]]
            obsInd2 = rev[rev[objInd2] : rev[objInd2 + 1]]
            a, b = esutil.numpy_util.match(data[self.visitKey][obsInd1], data[self.visitKey][obsInd2])
            matchedObsInd1.append(obsInd1[a])
            matchedObsInd2.append(obsInd2[b])
            matchedPairInd.append(np.full(len(a), ind))

        matchedObsInd1 = np.concatenate(matchedObsInd1)
        matchedObsInd2 = np.concatenate(matchedObsInd2)
        matchedPairInd = np.concatenate(matchedPairInd)

        separations = sphDist(
            np.deg2rad(np.array(data[self.raKey][matchedObsInd1])),
            np.deg2rad(np.array(data[self.decKey][matchedObsInd1])),
            np.deg2rad(np.array(data[self.raKey][matchedObsInd2])),
            np.deg2rad(np.array(data[self.decKey][matchedObsInd2])),
        )

        # Compute the mean from the ragged array of pairs by
        # using np.add.at to sum numerator and denominator.
        sepMean = np.zeros(len(i1))
        nSep = np.zeros_like(sepMean, dtype=np.int32)
        np.add.at(sepMean, matchedPairInd, separations)
        np.add.at(nSep, matchedPairInd, 1)
        good = nSep > 1
        sepMean[good] /= nSep[good]
        sepMean[~good] = np.nan

        # There are no good pairs, so return the default.
        if good.sum() == 0:
            return distanceParams

        # Compute the stdev with sqrt(sum((sep - mean(sep))**2.)/(nsep - 1))
        sepStd = np.zeros_like(sepMean)
        np.add.at(
            sepStd,
            matchedPairInd,
            (separations - sepMean[matchedPairInd]) ** 2.0,
        )
        sepStd[good] = np.sqrt(sepStd[good] / (nSep[good] - 1))
        rmsDistances = sepStd[good]

        # Need sepResiduals, but only when nSep is > 2.
        bad2 = nSep <= 2
        sepMean[bad2] = np.nan
        sepResiduals = separations - sepMean[matchedPairInd]
        sepResiduals = sepResiduals[np.isfinite(sepResiduals)]

        # This is always going to be valid because we checked the number
        # of good pairs above.
        AMx = (np.median(rmsDistances) * u.radian).to(u.marcsec)

        # Because there is a more stringent selection for sepResiduals,
        # we need to check that we have enough to compute the metrics.
        if len(sepResiduals) <= 1:
            AFx = np.nan * u.percent
            ADx = np.nan * u.marcsec
            absDiffSeparations = np.array([]) * u.marcsec
        else:
            absDiffSeparations = (abs(sepResiduals - np.median(sepResiduals)) * u.radian).to(u.marcsec)
            afThreshhold = 100.0 - self.threshAF
            ADx = np.percentile(absDiffSeparations, afThreshhold)
            AFx = 100 * np.mean(np.abs(absDiffSeparations) > self.threshAD * u.marcsec) * u.percent

        distanceParams["rmsDistances"] = (rmsDistances * u.radian).to(u.marcsec).value
        distanceParams["separationResiduals"] = absDiffSeparations.value
        distanceParams["AMx"] = AMx.value
        distanceParams["ADx"] = ADx.value
        distanceParams["AFx"] = AFx.value

        return distanceParams


def sphDist(ra_mean, dec_mean, ra, dec):
    """Calculate distance on the surface of a unit sphere.

    This function was borrowed from faro.

    Parameters
    ----------
    ra_mean : `float`
        Mean RA in radians.
    dec_mean : `float`
        Mean Dec in radians.
    ra : `numpy.array` [`float`]
        Array of RA in radians.
    dec : `numpy.array` [`float`]
        Array of Dec in radians.

    Notes
    -----
    Uses the Haversine formula to preserve accuracy at small angles.
    Law of cosines approach doesn't work well for the typically very small
    differences that we're looking at here.
    """
    # Haversine
    dra = ra - ra_mean
    ddec = dec - dec_mean
    a = np.square(np.sin(ddec / 2)) + np.cos(dec_mean) * np.cos(dec) * np.square(np.sin(dra / 2))
    dist = 2 * np.arcsin(np.sqrt(a))

    # This is what the law of cosines would look like
    #    dist = np.arccos(np.sin(dec1)*np.sin(dec2) +
    #                     np.cos(dec1)*np.cos(dec2)*np.cos(ra1 - ra2))

    # This will also work, but must run separately for each element
    # whereas the numpy version will run on either scalars or arrays:
    #   sp1 = geom.SpherePoint(ra1, dec1, geom.radians)
    #   sp2 = geom.SpherePoint(ra2, dec2, geom.radians)
    #   return sp1.separation(sp2).asRadians()

    return dist
