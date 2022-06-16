from __future__ import annotations

__all__ = (
    "StellarLocusFitAction"
)

import numpy as np
from scipy.stats import median_absolute_deviation as sigmaMad

from lsst.pex.config import DictField
from ..plotActions.plotUtils import stellarLocusFit, perpDistance
from ..interfaces import KeyedData, KeyedDataSchema, KeyedDataAction, Vector


class StellarLocusFitAction(KeyedDataAction):
    stellarLocusFitDict = DictField(
        doc="The parameters to use for the stellar locus fit. The default parameters are examples and are "
            "not useful for any of the fits. The dict needs to contain xMin/xMax/yMin/yMax which are the "
            "limits of the initial box for fitting the stellar locus, mHW and bHW are the initial "
            "intercept and gradient for the fitting.",
        keytype=str,
        itemtype=float,
        default={"xMin": 0.1, "xMax": 0.2, "yMin": 0.1, "yMax": 0.2, "mHW": 0.5, "bHW": 0.0}
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return (("x", Vector), ("y", Vector))

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result: KeyedData = {}

        xs = data['x']
        ys = data['y']
        fitParams = stellarLocusFit(xs, ys, self.stellarLocusFitDict)
        fitPoints = np.where((xs > fitParams["xMin"]) & (xs < fitParams["xMax"])
                             & (ys > fitParams["yMin"]) & (ys < fitParams["yMax"]))[0]

        if np.fabs(fitParams["mHW"]) > 1:
            ysFitLineHW = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLineHW = (ysFitLineHW - fitParams["bHW"])/fitParams["mHW"]
            ysFitLine = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLine = (ysFitLine - fitParams["bODR"])/fitParams["mODR"]
            ysFitLine2 = np.array([fitParams["yMin"], fitParams["yMax"]])
            xsFitLine2 = (ysFitLine2 - fitParams["bODR2"])/fitParams["mODR2"]

        else:
            xsFitLineHW = np.array([fitParams["xMin"], fitParams["xMax"]])
            ysFitLineHW = fitParams["mHW"]*xsFitLineHW + fitParams["bHW"]
            xsFitLine = [fitParams["xMin"], fitParams["xMax"]]
            ysFitLine = [fitParams["mODR"]*xsFitLine[0] + fitParams["bODR"],
                         fitParams["mODR"]*xsFitLine[1] + fitParams["bODR"]]
            xsFitLine2 = [fitParams["xMin"], fitParams["xMax"]]
            ysFitLine2 = [fitParams["mODR2"]*xsFitLine2[0] + fitParams["bODR2"],
                          fitParams["mODR2"]*xsFitLine2[1] + fitParams["bODR2"]]

        # Calculate the distances to that line
        # Need two points to characterise the lines we want
        # to get the distances to
        p1 = np.array([xsFitLine[0], ysFitLine[0]])
        p2 = np.array([xsFitLine[1], ysFitLine[1]])

        p1HW = np.array([xsFitLine[0], ysFitLineHW[0]])
        p2HW = np.array([xsFitLine[1], ysFitLineHW[1]])

        distsHW = perpDistance(p1HW, p2HW, zip(xs[fitPoints], ys[fitPoints]))
        dists = perpDistance(p1, p2, zip(xs[fitPoints], ys[fitPoints]))

        # Now we have the information for the perpendicular line we
        # can use it to calculate the points at the ends of the
        # perpendicular lines that intersect at the box edges
        if np.fabs(fitParams["mHW"]) > 1:
            xMid = (fitParams["yMin"] - fitParams["bODR2"])/fitParams["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"]*xs + fitParams["bPerpMin"]
        else:
            xs = np.array([fitParams["xMin"] - 0.2, fitParams["xMin"], fitParams["xMin"] + 0.2])
            ys = xs*fitParams["mPerp"] + fitParams["bPerpMin"]

        if np.fabs(fitParams["mHW"]) > 1:
            xMid = (fitParams["yMax"] - fitParams["bODR2"])/fitParams["mODR2"]
            xs = np.array([xMid - 0.5, xMid, xMid + 0.5])
            ys = fitParams["mPerp"]*xs + fitParams["bPerpMax"]
        else:
            xs = np.array([fitParams["xMax"] - 0.2, fitParams["xMax"], fitParams["xMax"] + 0.2])
            ys = xs*fitParams["mPerp"] + fitParams["bPerpMax"]

        result["{identifier}_sigmaMAD".format(**kwargs)] = sigmaMad(dists)
        result["{identifier}_median".format(**kwargs)] = np.median(dists)
        result["{identifier}_hardwired_sigmaMAD".format(**kwargs)] = sigmaMad(distsHW)
        result["{identifier}_hardwired_median".format(**kwargs)] = np.median(distsHW)

        return result
