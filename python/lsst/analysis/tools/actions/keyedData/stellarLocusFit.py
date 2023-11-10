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


def stellarLocusFit(xs, ys, paramDict):
    """Make a fit to the stellar locus.

    Parameters
    ----------
    xs : `numpy.ndarray`
        The color on the xaxis
    ys : `numpy.ndarray`
        The color on the yaxis
    paramDict : lsst.pex.config.dictField.Dict
        A dictionary of parameters for line fitting
        xMin : `float`
            The minimum x edge of the box to use for initial fitting
        xMax : `float`
            The maximum x edge of the box to use for initial fitting
        yMin : `float`
            The minimum y edge of the box to use for initial fitting
        yMax : `float`
            The maximum y edge of the box to use for initial fitting
        mHW : `float`
            The hardwired gradient for the fit
        bHW : `float`
            The hardwired intercept of the fit

    Returns
    -------
    paramsOut : `dict`
        A dictionary of the calculated fit parameters
        xMin : `float`
            The minimum x edge of the box to use for initial fitting
        xMax : `float`
            The maximum x edge of the box to use for initial fitting
        yMin : `float`
            The minimum y edge of the box to use for initial fitting
        yMax : `float`
            The maximum y edge of the box to use for initial fitting
        mHW : `float`
            The hardwired gradient for the fit
        bHW : `float`
            The hardwired intercept of the fit
        mODR : `float`
            The gradient calculated by the ODR fit
        bODR : `float`
            The intercept calculated by the ODR fit
        yBoxMin : `float`
            The y value of the fitted line at xMin
        yBoxMax : `float`
            The y value of the fitted line at xMax
        bPerpMin : `float`
            The intercept of the perpendicular line that goes through xMin
        bPerpMax : `float`
            The intercept of the perpendicular line that goes through xMax
        mODR2 : `float`
            The gradient from the second round of fitting
        bODR2 : `float`
            The intercept from the second round of fitting
        mPerp : `float`
            The gradient of the line perpendicular to the line from the
            second fit

    Notes
    -----
    The code does two rounds of fitting, the first is initiated using the
    hardwired values given in the `paramDict` parameter and is done using
    an Orthogonal Distance Regression fit to the points defined by the
    box of xMin, xMax, yMin and yMax. Once this fitting has been done a
    perpendicular bisector is calculated at either end of the line and
    only points that fall within these lines are used to recalculate the fit.
    """
    # Points to use for the fit
    fitPoints = np.where(
        (xs > paramDict["xMin"])
        & (xs < paramDict["xMax"])
        & (ys > paramDict["yMin"])
        & (ys < paramDict["yMax"])
    )[0]

    linear = scipyODR.polynomial(1)

    data = scipyODR.Data(xs[fitPoints], ys[fitPoints])
    odr = scipyODR.ODR(data, linear, beta0=[paramDict["bHW"], paramDict["mHW"]])
    params = odr.run()
    mODR = float(params.beta[1])
    bODR = float(params.beta[0])

    paramsOut = {
        "xMin": paramDict["xMin"],
        "xMax": paramDict["xMax"],
        "yMin": paramDict["yMin"],
        "yMax": paramDict["yMax"],
        "mHW": paramDict["mHW"],
        "bHW": paramDict["bHW"],
        "mODR": mODR,
        "bODR": bODR,
    }

    # Having found the initial fit calculate perpendicular ends
    mPerp = -1.0 / mODR
    # When the gradient is really steep we need to use
    # the y limits of the box rather than the x ones

    if np.abs(mODR) > 1:
        yBoxMin = paramDict["yMin"]
        xBoxMin = (yBoxMin - bODR) / mODR
        yBoxMax = paramDict["yMax"]
        xBoxMax = (yBoxMax - bODR) / mODR
    else:
        yBoxMin = mODR * paramDict["xMin"] + bODR
        xBoxMin = paramDict["xMin"]
        yBoxMax = mODR * paramDict["xMax"] + bODR
        xBoxMax = paramDict["xMax"]

    bPerpMin = yBoxMin - mPerp * xBoxMin

    paramsOut["yBoxMin"] = yBoxMin
    paramsOut["bPerpMin"] = bPerpMin

    bPerpMax = yBoxMax - mPerp * xBoxMax

    paramsOut["yBoxMax"] = yBoxMax
    paramsOut["bPerpMax"] = bPerpMax

    # Use these perpendicular lines to chose the data and refit
    fitPoints = (ys > mPerp * xs + bPerpMin) & (ys < mPerp * xs + bPerpMax)
    data = scipyODR.Data(xs[fitPoints], ys[fitPoints])
    odr = scipyODR.ODR(data, linear, beta0=[bODR, mODR])
    params = odr.run()
    mODR = float(params.beta[1])
    bODR = float(params.beta[0])

    paramsOut["mODR2"] = float(params.beta[1])
    paramsOut["bODR2"] = float(params.beta[0])

    paramsOut["mPerp"] = -1.0 / paramsOut["mODR2"]

    return paramsOut


def perpDistance(p1, p2, points):
    """Calculate the perpendicular distance to a line from a point.

    Parameters
    ----------
    p1 : `numpy.ndarray`
        A point on the line
    p2 : `numpy.ndarray`
        Another point on the line
    points : `zip`
        The points to calculate the distance to

    Returns
    -------
    dists : `list`
        The distances from the line to the points. Uses the cross
        product to work this out.
    """
    dists = []
    for point in points:
        point = np.array(point)
        distToLine = np.cross(p2 - p1, point - p1) / np.linalg.norm(p2 - p1)
        dists.append(distToLine)

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
            (f"{self.identity or ''}_hardwired_sigmaMAD", Scalar),
            (f"{self.identity or ''}_hardwired_median", Scalar),
        )
        return value

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        xs = cast(Vector, data["x"])
        ys = cast(Vector, data["y"])
        fitParams = stellarLocusFit(xs, ys, self.stellarLocusFitDict)
        fitPoints = np.where(
            (xs > fitParams["xMin"])  # type: ignore
            & (xs < fitParams["xMax"])  # type: ignore
            & (ys > fitParams["yMin"])  # type: ignore
            & (ys < fitParams["yMax"])  # type: ignore
        )[0]

        if np.fabs(fitParams["mHW"]) > 1:
            ysFitLineHW = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLineHW = (ysFitLineHW - fitParams["bHW"]) / fitParams["mHW"]
            ysFitLine = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLine = (ysFitLine - fitParams["bODR"]) / fitParams["mODR"]
            ysFitLine2 = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLine2 = (ysFitLine2 - fitParams["bODR2"]) / fitParams["mODR2"]

        else:
            xsFitLineHW = np.array([fitParams["xMin"], fitParams["xMax"]])
            ysFitLineHW = fitParams["mHW"] * xsFitLineHW + fitParams["bHW"]
            xsFitLine = [fitParams["xMin"], fitParams["xMax"]]
            ysFitLine = np.array(
                [
                    fitParams["mODR"] * xsFitLine[0] + fitParams["bODR"],
                    fitParams["mODR"] * xsFitLine[1] + fitParams["bODR"],
                ]
            )
            xsFitLine2 = [fitParams["xMin"], fitParams["xMax"]]
            ysFitLine2 = np.array(
                [
                    fitParams["mODR2"] * xsFitLine2[0] + fitParams["bODR2"],
                    fitParams["mODR2"] * xsFitLine2[1] + fitParams["bODR2"],
                ]
            )

        # Calculate the distances to that line
        # Need two points to characterise the lines we want
        # to get the distances to
        p1 = np.array([xsFitLine[0], ysFitLine[0]])
        p2 = np.array([xsFitLine[1], ysFitLine[1]])

        p1HW = np.array([xsFitLine[0], ysFitLineHW[0]])
        p2HW = np.array([xsFitLine[1], ysFitLineHW[1]])

        # Convert this to mmag
        distsHW = np.array(perpDistance(p1HW, p2HW, zip(xs[fitPoints], ys[fitPoints]))) * 1000
        dists = np.array(perpDistance(p1, p2, zip(xs[fitPoints], ys[fitPoints]))) * 1000

        # Now we have the information for the perpendicular line we
        # can use it to calculate the points at the ends of the
        # perpendicular lines that intersect at the box edges
        if np.fabs(fitParams["mHW"]) > 1:
            xMid = (fitParams["yMin"] - fitParams["bODR2"]) / fitParams["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"] * xs + fitParams["bPerpMin"]
        else:
            xs = np.array([fitParams["xMin"] - 0.2, fitParams["xMin"], fitParams["xMin"] + 0.2])
            ys = xs * fitParams["mPerp"] + fitParams["bPerpMin"]

        if np.fabs(fitParams["mHW"]) > 1:
            xMid = (fitParams["yMax"] - fitParams["bODR2"]) / fitParams["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"] * xs + fitParams["bPerpMax"]
        else:
            xs = np.array([fitParams["xMax"] - 0.2, fitParams["xMax"], fitParams["xMax"] + 0.2])
            ys = xs * fitParams["mPerp"] + fitParams["bPerpMax"]

        fit_sigma, fit_med = (sigmaMad(dists), np.median(dists)) if len(dists) else (np.nan, np.nan)
        fitParams[f"{self.identity or ''}_sigmaMAD"] = fit_sigma
        fitParams[f"{self.identity or ''}_median"] = fit_med
        fit_sigma, fit_med = (sigmaMad(distsHW), np.median(distsHW)) if len(distsHW) else (np.nan, np.nan)
        fitParams[f"{self.identity or ''}_hardwired_sigmaMAD"] = fit_sigma
        fitParams[f"{self.identity or ''}_hardwired_median"] = fit_med

        return fitParams  # type: ignore
