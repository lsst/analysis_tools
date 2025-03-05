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

__all__ = (
    "getTractCorners",
    "getPatchCorners",
    "http_client",
)

from collections.abc import Generator
from contextlib import contextmanager

import numpy as np
import requests
from lsst.geom import Box2D
from requests.adapters import HTTPAdapter
from urllib3 import Retry


def getTractCorners(skymap, tractId):
    """Calculate the corners of a tract, given a skymap.

    Parameters
    ----------
    skymap : `lsst.skymap`
    tractId : `int`
        Identification number of the tract whose corner coordinates
        are returned.

    Returns
    -------
    corners : `list` of `tuples` of `float`

    Notes
    -----
    Corners are returned in degrees and wrapped in ra.
    """
    tractCorners = skymap[tractId].getVertexList()
    corners = _wrapRa([(corner.getRa().asDegrees(), corner.getDec().asDegrees()) for corner in tractCorners])

    return corners


def getPatchCorners(tractInfo, patchId):
    """Calculate the corners of a patch, given tractInfo.

    Parameters
    ----------
    tractInfo : `lsst.skymap.tractInfo.ExplicitTractInfo`
        Tract info object of the tract containing the patch whose
        corner coordinates are returned.
    patchId : `int`
        Identification number of the patch whose corner coordinates
        are returned.

    Returns
    -------
    corners : `list` of `tuples` of `float`

    Notes
    -----
    Corners are returned in degrees and are wrapped in ra.
    """
    patchInfo = tractInfo.getPatchInfo(patchId)
    patchCorners = Box2D(patchInfo.getInnerBBox()).getCorners()

    tractWcs = tractInfo.getWcs()
    patchCorners = tractWcs.pixelToSky(patchCorners)
    corners = _wrapRa([(corner.getRa().asDegrees(), corner.getDec().asDegrees()) for corner in patchCorners])

    return corners


def _wrapRa(corners):
    """Wrap in right ascension if the corners span RA=0

    Parameters
    ----------
    corners :  `list` of `tuples` of `float`
        Pairs of coordinates representing tract or patch corners.

    Returns
    -------
    corners : `list` of `tuples` of `float`
        Pairs of coordinates representing tract or patch corners,
        wrapped in RA.
    """

    minRa = np.min([corner[0] for corner in corners])
    maxRa = np.max([corner[0] for corner in corners])
    # If the tract needs wrapping in ra, wrap it
    if maxRa - minRa > 10:
        x = maxRa
        maxRa = 360 + minRa
        minRa = x
        minDec = np.min([corner[1] for corner in corners])
        maxDec = np.max([corner[1] for corner in corners])
        corners = [(minRa, minDec), (maxRa, minDec), (maxRa, maxDec), (minRa, maxDec)]

    return corners


@contextmanager
def http_client() -> Generator[requests.Session]:
    """Creates a requests session with a custom transport to support
    automatic retries with backoff for dealing with transient server-side
    issues.

    Notes
    -----
    The goal of the adapter defined here is to avoid premature client abends
    when transient server or infrastructure issues prevent good-faith attempts
    at accessing APIs. To the extent that we want to balance "eventually
    successful" HTTP requests with the desire to vacate the compute resources
    our process is occupying, these retries should not overstay their welcome.

    The "POST" HTTP verb is not usually part of the allowed methods for retries
    because unlike "PUT", "POST" is not considered idempotent by default. It is
    partially for this reason that a custom Retry adapter is needed, because
    by default "POST" requests would not be retried for status.

    The backoff_factor is an exponential factor used to calculate how long to
    sleep between the third and subsequent tries, in seconds. The first retry
    is immediate and the total backoff won't exceed backoff_max, which defaults
    to 120 seconds.
    """

    retriable_statuses = [
        requests.codes.too_many_requests,
        requests.codes.server_error,
        requests.codes.bad_gateway,
        requests.codes.service_unavailable,
        requests.codes.gateway_timeout,
    ]
    session = requests.Session()
    retry_strategy = Retry(
        total=None,  # use specific conditional constraints
        connect=3,  # network or tcp errors
        read=0,  # request sent, response is bad
        status=5,  # retries based on bad response status (see retriable_statuses)
        redirect=3,  # default value, follow 3 redirects
        other=0,  # edge cases and weird stuff
        backoff_factor=0.1,  # sleep == {factor} * 2^(previous tries)
        status_forcelist=retriable_statuses,
        raise_on_status=True,
        allowed_methods={"GET", "HEAD", "POST", "PUT"},
    )
    session.mount("http://", HTTPAdapter(max_retries=retry_strategy))
    session.mount("https://", HTTPAdapter(max_retries=retry_strategy))
    try:
        yield session
    finally:
        session.close()
