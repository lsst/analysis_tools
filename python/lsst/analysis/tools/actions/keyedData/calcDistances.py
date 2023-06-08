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
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
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
        data: KeyedData
            Catalog of data including coordinate, visit, and object group
            information.
        Returns
        -------
        distanceParams: `dict`
            Dictionary of the calculated arrays and metrics with the following
            keys:
            - ``rmsDistances`` : Per-object rms of separations (`np.array`).
            - ``separationResiduals`` : All separations minus per-object median
                (`np.array`)
            - ``AMx`` : AMx metric (`float`).
            - ``ADx`` : ADx metric (`float`).
            - ``AFx`` : AFx metric (`float`).
        """
        D = self.annulus * u.arcmin
        width = self.width * u.arcmin
        annulus = (D + (width / 2) * np.array([-1, +1])).to(u.radian)

        df = pd.DataFrame(
            {
                "groupKey": data[self.groupKey],
                "coord_ra": data[self.raKey],
                "coord_dec": data[self.decKey],
                "visit": data[self.visitKey],
            }
        )

        meanRa = df.groupby("groupKey")["coord_ra"].aggregate("mean")
        meanDec = df.groupby("groupKey")["coord_dec"].aggregate("mean")

        catalog = SkyCoord(meanRa.to_numpy() * u.degree, meanDec.to_numpy() * u.degree)
        idx, idxCatalog, d2d, d3d = catalog.search_around_sky(catalog, annulus[1])
        inAnnulus = d2d > annulus[0]
        idx = idx[inAnnulus]
        idxCatalog = idxCatalog[inAnnulus]

        rmsDistances = []
        sepResiduals = []
        for id in range(len(meanRa)):
            match_inds = idx == id
            match_ids = idxCatalog[match_inds & (idxCatalog != id)]
            if match_ids.sum() == 0:
                continue

            object_srcs = df.loc[df["groupKey"] == meanRa.index[id]]

            object_visits = object_srcs["visit"].to_numpy()
            object_ras = (object_srcs["coord_ra"].to_numpy() * u.degree).to(u.radian).value
            object_decs = (object_srcs["coord_dec"].to_numpy() * u.degree).to(u.radian).value
            if len(object_srcs) <= 1:
                continue
            object_srcs = object_srcs.set_index("visit")
            object_srcs.sort_index(inplace=True)

            for id2 in match_ids:
                match_srcs = df.loc[df["groupKey"] == meanRa.index[id2]]
                match_visits = match_srcs["visit"].to_numpy()
                match_ras = (match_srcs["coord_ra"].to_numpy() * u.degree).to(u.radian).value
                match_decs = (match_srcs["coord_dec"].to_numpy() * u.degree).to(u.radian).value

                separations = matchVisitComputeDistance(
                    object_visits, object_ras, object_decs, match_visits, match_ras, match_decs
                )

                if len(separations) > 1:
                    rmsDist = np.std(separations, ddof=1)
                    rmsDistances.append(rmsDist)
                if len(separations) > 2:
                    sepResiduals.append(separations - np.median(separations))

        if len(rmsDistances) == 0:
            AMx = np.nan * u.marcsec
        else:
            AMx = (np.median(rmsDistances) * u.radian).to(u.marcsec)

        if len(sepResiduals) <= 1:
            AFx = np.nan * u.percent
            ADx = np.nan * u.marcsec
            absDiffSeparations = np.array([]) * u.marcsec
        else:
            sepResiduals = np.concatenate(sepResiduals)
            absDiffSeparations = (abs(sepResiduals - np.median(sepResiduals)) * u.radian).to(u.marcsec)
            afThreshhold = 100.0 - self.threshAF
            ADx = np.percentile(absDiffSeparations, afThreshhold)
            AFx = 100 * np.mean(np.abs(absDiffSeparations) > self.threshAD * u.marcsec) * u.percent

        distanceParams = {
            "rmsDistances": (np.array(rmsDistances) * u.radian).to(u.marcsec).value,
            "separationResiduals": absDiffSeparations.value,
            "AMx": AMx.value,
            "ADx": ADx.value,
            "AFx": AFx.value,
        }

        return distanceParams


def matchVisitComputeDistance(visit_obj1, ra_obj1, dec_obj1, visit_obj2, ra_obj2, dec_obj2):
    """Calculate obj1-obj2 distance for each visit in which both objects are
    seen.

    For each visit shared between visit_obj1 and visit_obj2, calculate the
    spherical distance between the obj1 and obj2. visit_obj1 and visit_obj2 are
    assumed to be unsorted. This function was borrowed from faro.

    Parameters
    ----------
    visit_obj1 : scalar, list, or numpy.array of int or str
        List of visits for object 1.
    ra_obj1 : scalar, list, or numpy.array of float
        List of RA in each visit for object 1.  [radians]
    dec_obj1 : scalar, list or numpy.array of float
        List of Dec in each visit for object 1. [radians]
    visit_obj2 : list or numpy.array of int or str
        List of visits for object 2.
    ra_obj2 : list or numpy.array of float
        List of RA in each visit for object 2.  [radians]
    dec_obj2 : list or numpy.array of float
        List of Dec in each visit for object 2.  [radians]
    Results
    -------
    list of float
        spherical distances (in radians) for matching visits.
    """
    distances = []
    visit_obj1_idx = np.argsort(visit_obj1)
    visit_obj2_idx = np.argsort(visit_obj2)
    j_raw = 0
    j = visit_obj2_idx[j_raw]
    for i in visit_obj1_idx:
        while (visit_obj2[j] < visit_obj1[i]) and (j_raw < len(visit_obj2_idx) - 1):
            j_raw += 1
            j = visit_obj2_idx[j_raw]
        if visit_obj2[j] == visit_obj1[i]:
            if np.isfinite([ra_obj1[i], dec_obj1[i], ra_obj2[j], dec_obj2[j]]).all():
                distances.append(sphDist(ra_obj1[i], dec_obj1[i], ra_obj2[j], dec_obj2[j]))
    return distances


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
