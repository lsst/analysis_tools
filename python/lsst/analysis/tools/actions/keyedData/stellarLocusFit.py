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

__all__ = ("StellarLocusFitAction",)

from typing import cast

import numpy as np
import scipy.odr as scipyODR
from lsst.pex.config import DictField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Scalar, Vector
from ...math import sigmaMad


def stellarLocusFit(xs, ys, mags, paramDict):
    """Make a fit to the stellar locus.

    Parameters
    ----------
    xs : `numpy.ndarray` [`float`]
        The color on the xaxis.
    ys : `numpy.ndarray` [`float`]
        The color on the yaxis.
    mags : `numpy.ndarray` [`float`]
        The magnitude of the reference band flux (in mag).
    paramDict : `dict` [`str`, `float`]
        A dictionary of parameters for line fitting:

        ``"xMin"``
            The minimum x edge of the box to use for initial fitting (`float`).
        ``"xMax"``
            The maximum x edge of the box to use for initial fitting (`float`).
        ``"yMin"``
            The minimum y edge of the box to use for initial fitting (`float`).
        ``"yMax"``
            The maximum y edge of the box to use for initial fitting (`float`).
        ``"mHW"``
            The hardwired gradient for the fit (`float`).
        ``"bHw"``
            The hardwired intercept of the fit (`float`).

    Returns
    -------
    fitParams : `dict`
        A dictionary of the calculated fit parameters:

        ``"bPerpMin"``
            The intercept of the perpendicular line that goes through xMin
            (`float`).
        ``"bPerpMax"``
            The intercept of the perpendicular line that goes through xMax
            (`float`).
        ``"mODR"``
            The gradient from the final round of fitting (`float`).
        ``"bODR"``
            The intercept from the final round of fitting (`float`).
        ``"mPerp"``
            The gradient of the line perpendicular to the line from the final
            fit (`float`).
        ``"fitPoints"``
            A boolean list indicating which points were used in the final fit
            (`list` [`bool`]).

    Notes
    -----
    The code does two rounds of fitting, the first is initiated using the
    hardwired values given in ``paramDict`` and is done using an Orthogonal
    Distance Regression (ODR) fit to the points defined by the box with limits
    defined by the keys: xMin, xMax, yMin, and yMax. Once this fitting has been
    done a perpendicular bisector is calculated at either end of the line and
    only points that fall within these lines are used to recalculate the fit.
    We also perform clipping of points perpendicular to the fit line that have
    distances that deviate more than 5-sigma from the fit.
    """
    fitParams = {}
    # Initial subselection of points to use for the fit
    # Check for nans/infs
    goodPoints = np.isfinite(xs) & np.isfinite(ys) & np.isfinite(mags)

    fitPoints = (
        goodPoints
        & (xs > paramDict["xMin"])
        & (xs < paramDict["xMax"])
        & (ys > paramDict["yMin"])
        & (ys < paramDict["yMax"])
    )
    linear = scipyODR.polynomial(1)

    fitData = scipyODR.Data(xs[fitPoints], ys[fitPoints])
    odr = scipyODR.ODR(fitData, linear, beta0=[paramDict["bHW"], paramDict["mHW"]])
    params = odr.run()
    mODR0 = float(params.beta[1])
    bODR0 = float(params.beta[0])

    # Having found the initial fit calculate perpendicular ends
    mPerp0 = -1.0 / mODR0

    # When the gradient is really steep we need to use
    # the y limits of the fit line rather than the x ones.

    if np.abs(mODR0) > 1:
        yPerpMin = paramDict["yMin"]
        xPerpMin = (yPerpMin - bODR0) / mODR0
        yPerpMax = paramDict["yMax"]
        xPerpMax = (yPerpMax - bODR0) / mODR0
    else:
        yPerpMin = mODR0 * paramDict["xMin"] + bODR0
        xPerpMin = paramDict["xMin"]
        yPerpMax = mODR0 * paramDict["xMax"] + bODR0
        xPerpMax = paramDict["xMax"]

    bPerpMin = yPerpMin - mPerp0 * xPerpMin
    bPerpMax = yPerpMax - mPerp0 * xPerpMax

    fitParams["bPerpMin"] = bPerpMin
    fitParams["bPerpMax"] = bPerpMax

    fitParams["mODR"] = float(params.beta[1])
    fitParams["bODR"] = float(params.beta[0])

    fitParams["mPerp"] = -1.0 / fitParams["mODR"]
    fitParams["goodPoints"] = goodPoints
    fitParams["fitPoints"] = fitPoints
    fitParams["paramDict"] = paramDict

    return fitParams


def perpDistance(p1, p2, points):
    """Calculate the perpendicular distance to a line from a point.

    Parameters
    ----------
    p1 : `numpy.ndarray` [`float`]
        A point on the line.
    p2 : `numpy.ndarray` [`float`]
        Another point on the line.
    points : `zip` [(`float`, `float`)]
        The points to calculate the distance to.

    Returns
    -------
    dists : `numpy.ndarray` [`float`]
        The distances from the line to the points. Uses the cross
        product to work this out.
    """
    if sum(p2 - p1) == 0:
        raise ValueError(f"Must supply two different points for p1, p2. Got {p1}, {p2}")
    points = list(points)
    if len(points) == 0:
        raise ValueError("Must provied a non-empty zip() list of points.")
    dists = np.cross(p2 - p1, points - p1) / np.linalg.norm(p2 - p1)

    return dists


class StellarLocusFitAction(KeyedDataAction):
    r"""Determine Stellar Locus fit parameters from given input `Vector`\ s."""

    stellarLocusFitDict = DictField[str, float](
        doc="The parameters to use for the stellar locus fit. The default parameters are examples and are "
        "not useful for any of the fits. The dict needs to contain xMin/xMax/yMin/yMax which are the "
        "limits of the initial box for fitting the stellar locus, mHW and bHW are the initial "
        "intercept and gradient for the fitting.",
        default={"xMin": 0.1, "xMax": 0.2, "yMin": 0.1, "yMax": 0.2, "mHW": 0.5, "bHW": 0.0},
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return (("x", Vector), ("y", Vector))

    def getOutputSchema(self) -> KeyedDataSchema:
        value = (
            (f"{self.identity or ''}_sigmaMAD", Scalar),
            (f"{self.identity or ''}_median", Scalar),
        )
        return value

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        xs = cast(Vector, data["x"])
        ys = cast(Vector, data["y"])
        mags = cast(Vector, data["mag"])

        fitParams = stellarLocusFit(xs, ys, mags, self.stellarLocusFitDict)
        fitPoints = fitParams["fitPoints"]

        if np.fabs(self.stellarLocusFitDict["mHW"]) > 1:
            ysFitLineHW = np.array([self.stellarLocusFitDict["yMin"], self.stellarLocusFitDict["yMax"]])
            xsFitLineHW = (ysFitLineHW - self.stellarLocusFitDict["bHW"]) / self.stellarLocusFitDict["mHW"]
            ysFitLine = np.array([self.stellarLocusFitDict["yMin"], self.stellarLocusFitDict["yMax"]])
            xsFitLine = (ysFitLine - fitParams["bODR"]) / fitParams["mODR"]

        else:
            xsFitLineHW = np.array([self.stellarLocusFitDict["xMin"], self.stellarLocusFitDict["xMax"]])
            ysFitLineHW = self.stellarLocusFitDict["mHW"] * xsFitLineHW + self.stellarLocusFitDict["bHW"]
            xsFitLine = [self.stellarLocusFitDict["xMin"], self.stellarLocusFitDict["xMax"]]
            ysFitLine = np.array(
                [
                    fitParams["mODR"] * xsFitLine[0] + fitParams["bODR"],
                    fitParams["mODR"] * xsFitLine[1] + fitParams["bODR"],
                ]
            )

        # Calculate the distances to that line.
        # Need two points to characterize the lines we want to get the
        # distances to.
        p1 = np.array([xsFitLine[0], ysFitLine[0]])
        p2 = np.array([xsFitLine[1], ysFitLine[1]])

        # Convert this to mmag.
        dists = np.array(perpDistance(p1, p2, zip(xs[fitPoints], ys[fitPoints]))) * 1000

        # Now we have the information for the perpendicular line we
        # can use it to calculate the points at the ends of the
        # perpendicular lines that intersect at the box edges.
        if np.fabs(self.stellarLocusFitDict["mHW"]) > 1:
            xMid = (self.stellarLocusFitDict["yMin"] - fitParams["bODR"]) / fitParams["mODR"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"] * xs + fitParams["bPerpMin"]
        else:
            xs = np.array(
                [
                    self.stellarLocusFitDict["xMin"] - 0.2,
                    self.stellarLocusFitDict["xMin"],
                    self.stellarLocusFitDict["xMin"] + 0.2,
                ]
            )
            ys = xs * fitParams["mPerp"] + fitParams["bPerpMin"]

        if np.fabs(self.stellarLocusFitDict["mHW"]) > 1:
            xMid = (self.stellarLocusFitDict["yMax"] - fitParams["bODR"]) / fitParams["mODR"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"] * xs + fitParams["bPerpMax"]
        else:
            xs = np.array(
                [
                    self.stellarLocusFitDict["xMax"] - 0.2,
                    self.stellarLocusFitDict["xMax"],
                    self.stellarLocusFitDict["xMax"] + 0.2,
                ]
            )
            ys = xs * fitParams["mPerp"] + fitParams["bPerpMax"]

        fit_sigma, fit_med = (sigmaMad(dists), np.median(dists)) if len(dists) else (np.nan, np.nan)
        fitParams[f"{self.identity or ''}_sigmaMAD"] = fit_sigma
        fitParams[f"{self.identity or ''}_median"] = fit_med

        return fitParams  # type: ignore
