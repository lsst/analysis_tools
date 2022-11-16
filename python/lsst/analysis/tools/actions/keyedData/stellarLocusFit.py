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
from lsst.pex.config import DictField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Scalar, Vector
from ...statistics import sigmaMad
from ..plot.plotUtils import perpDistance, stellarLocusFit


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

        fitParams[f"{self.identity or ''}_sigmaMAD"] = sigmaMad(dists)
        fitParams[f"{self.identity or ''}_median"] = np.median(dists)
        fitParams[f"{self.identity or ''}_hardwired_sigmaMAD"] = sigmaMad(distsHW)
        fitParams[f"{self.identity or ''}_hardwired_median"] = np.median(distsHW)

        return fitParams  # type: ignore
