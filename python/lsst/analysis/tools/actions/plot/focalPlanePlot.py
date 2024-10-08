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

__all__ = ("FocalPlanePlot", "FocalPlaneGeometryPlot")

from typing import Mapping, Optional

import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS, Camera
from lsst.pex.config import ChoiceField, Field
from matplotlib.collections import PatchCollection
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import binned_statistic_2d, binned_statistic_dd

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from ...math import nanMax, nanMedian, nanMin, nanSigmaMad
from .plotUtils import addPlotInfo, mkColormap, sortAllArrays


class FocalPlanePlot(PlotAction):
    """Plots the focal plane distribution of a parameter.

    Given the detector positions in x and y, the focal plane positions are
    calculated using the camera model. A 2d binned statistic (default is mean)
    is then calculated and plotted for the parameter z as a function of the
    focal plane coordinates.
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis.", default="x (mm)", optional=True)
    yAxisLabel = Field[str](doc="Label to use for the y axis.", default="y (mm)", optional=True)
    zAxisLabel = Field[str](doc="Label to use for the z axis.", optional=False)
    nBins = Field[int](
        doc="Number of bins to use within the effective plot ranges along the spatial directions.",
        default=200,
    )
    doUseAdaptiveBinning = Field[bool](
        doc="If set to True, the number of bins is adapted to the source"
        " density, with lower densities using fewer bins. Under these"
        " circumstances the nBins parameter sets the minimum number of bins.",
        default=False,
    )
    statistic = Field[str](
        doc="Operation to perform in binned_statistic_2d",
        default="mean",
    )
    plotMin = Field[float](
        doc="Minimum in z-value to display in the focal plane plot and in the histogram plot, if applicable",
        default=None,
        optional=True,
    )
    plotMax = Field[float](
        doc="Maximum in z-value to display in the focal plane plot and in the histogram plot, if applicable",
        default=None,
        optional=True,
    )
    showStats = Field[bool](doc="Show statistics for plotted data", default=True)
    addHistogram = Field[bool](doc="Add a histogram of all input points", default=False)
    histBins = Field[int](doc="Number of bins to use in histogram", default=30)

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)
        # table is a dict that needs: x, y, run, skymap, filter, tract,

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """NOTE currently can only check that something is not a Scalar, not
        check that the data is consistent with Vector
        """
        needed = self.getInputSchema(**kwargs)
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in input")
        for name, typ in needed:
            isScalar = issubclass((colType := type(data[name.format(**kwargs)])), Scalar)
            if isScalar and typ != Scalar:
                raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}")

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        base.append(("x", Vector))
        base.append(("y", Vector))
        base.append(("z", Vector))
        base.append(("statMask", Vector))

        return base

    def statsAndText(self, arr, mask=None):
        """Calculate some stats from an array and return them
        and some text.
        """
        numPoints = len(arr)
        if mask is not None:
            arr = arr[mask]
        med = nanMedian(arr)
        sigMad = nanSigmaMad(arr)

        statsText = (
            "Median: {:0.2f}\n".format(med)
            + r"$\sigma_{MAD}$: "
            + "{:0.2f}\n".format(sigMad)
            + r"n$_{points}$: "
            + "{}".format(numPoints)
        )

        return med, sigMad, statsText

    def _addHistogram(self, histAx, data):
        bins = np.linspace(
            (self.plotMin if self.plotMin else nanMin(data.astype(float))),
            (self.plotMax if self.plotMax else nanMax(data.astype(float))),
            self.histBins,
        )
        histAx.hist(data.astype(float), bins=bins)
        histAx.set_xlabel(self.zAxisLabel)
        histAx.set_ylabel("Bin count")
        underflow = np.count_nonzero(data < bins[0])
        overflow = np.count_nonzero(data > bins[-1])
        nonfinite = np.count_nonzero(~np.isfinite(data))
        text = f"Underflow = {underflow}\nOverflow = {overflow}\nNon-Finite = {nonfinite}"
        anchored_text = AnchoredText(text, loc=1, pad=0.5)
        histAx.add_artist(anchored_text)

    def makePlot(
        self,
        data: KeyedData,
        camera: Camera,
        plotInfo: Optional[Mapping[str, str]] = None,
        **kwargs,
    ) -> Figure:
        """Prep the catalogue and then make a focalPlanePlot of the given
        column.

        Uses the axisLabels config options `x` and `y` to make an image, where
        the color corresponds to the 2d binned statistic (the mean is the
        default) applied to the `z` column. A summary panel is shown in the
        upper right corner of the resultant plot. The code uses the
        selectorActions to decide which points to plot and the
        statisticSelector actions to determine which points to use for the
        printed statistics.

        Parameters
        ----------
        data : `pandas.core.frame.DataFrame`
            The catalog to plot the points from.
        camera : `lsst.afw.cameraGeom.Camera`
            The camera used to map from pixel to focal plane positions.
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:

            ``"run"``
                The output run for the plots (`str`).
            ``"skymap"``
                The type of skymap used for the data (`str`).
            ``"filter"``
                The filter used for this data (`str`).
            ``"tract"``
                The tract that the data comes from (`str`).
            ``"bands"``
                The band(s) that the data comes from (`list` of `str`).

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.
        """
        if plotInfo is None:
            plotInfo = {}

        # This is a little hacky but is the easiest way to solve
        # the bands being printed on the plot correctly.
        if "band" in plotInfo.keys():
            plotInfo["bands"] = [plotInfo["band"]]

        if len(data["x"]) == 0:
            noDataFig = Figure()
            noDataFig.text(0.3, 0.5, "No data to plot after selectors applied")
            noDataFig = addPlotInfo(noDataFig, plotInfo)
            return noDataFig

        if self.addHistogram:
            fig, [ax, histAx] = plt.subplots(1, 2, dpi=300, figsize=(12, 6), width_ratios=[3, 2])
        else:
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(111)

        detectorIds = np.unique(data["detector"])
        focalPlane_x = np.zeros(len(data["x"]))
        focalPlane_y = np.zeros(len(data["y"]))
        for detectorId in detectorIds:
            detector = camera[detectorId]
            map = detector.getTransform(PIXELS, FOCAL_PLANE).getMapping()

            detectorInd = data["detector"] == detectorId
            points = np.array([data["x"][detectorInd], data["y"][detectorInd]])

            fp_x, fp_y = map.applyForward(points)
            focalPlane_x[detectorInd] = fp_x
            focalPlane_y[detectorInd] = fp_y

        if self.doUseAdaptiveBinning:
            # Use a course 32x32 binning to determine the mean source density
            # in regions where there are sources.
            binsx = np.linspace(focalPlane_x.min() - 1e-5, focalPlane_x.max() + 1e-5, 33)
            binsy = np.linspace(focalPlane_y.min() - 1e-5, focalPlane_y.max() + 1e-5, 33)

            binnedNumSrc = np.histogram2d(focalPlane_x, focalPlane_y, bins=[binsx, binsy])[0]
            meanSrcDensity = np.mean(binnedNumSrc, where=binnedNumSrc > 0.0)

            numBins = int(np.round(16.0 * np.sqrt(meanSrcDensity)))
            numBins = max(numBins, self.nBins)
        else:
            numBins = self.nBins

        # Add an arbitrary small offset to bins to ensure that the minimum does
        # not equal the maximum.
        binsx = np.linspace(focalPlane_x.min() - 1e-5, focalPlane_x.max() + 1e-5, numBins)
        binsy = np.linspace(focalPlane_y.min() - 1e-5, focalPlane_y.max() + 1e-5, numBins)

        statistic, x_edge, y_edge, binnumber = binned_statistic_2d(
            focalPlane_x, focalPlane_y, data["z"], statistic=self.statistic, bins=[binsx, binsy]
        )
        binExtent = [x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]]

        if self.showStats:
            sortedArrs = sortAllArrays([data["z"], data["x"], data["y"], data["statMask"]])
            [colorVals, xs, ys, stat] = sortedArrs
            statMed, statMad, statsText = self.statsAndText(colorVals, mask=stat)
            bbox = dict(facecolor="paleturquoise", alpha=0.5, edgecolor="none")
            ax.text(0.8, 0.91, statsText, transform=fig.transFigure, fontsize=8, bbox=bbox)

        median = nanMedian(statistic.ravel())
        mad = nanSigmaMad(statistic.ravel())

        vmin = self.plotMin if (self.plotMin is not None) else (median - 2 * mad)
        vmax = self.plotMax if (self.plotMax is not None) else (median + 2 * mad)

        plot = ax.imshow(statistic.T, extent=binExtent, vmin=vmin, vmax=vmax, origin="lower")

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(plot, cax=cax, extend="both")
        text = cax.text(
            0.5,
            0.5,
            self.zAxisLabel,
            color="k",
            rotation="vertical",
            transform=cax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )
        text.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])
        cax.tick_params(labelsize=7)

        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)
        ax.tick_params(axis="x", labelrotation=25)
        ax.tick_params(labelsize=7)

        ax.set_aspect("equal")

        if self.addHistogram:
            self._addHistogram(histAx, data["z"])

        plt.draw()

        # Add useful information to the plot
        plt.subplots_adjust(left=0.05, right=0.95)
        fig = plt.gcf()
        if plotInfo:
            fig = addPlotInfo(fig, plotInfo)

        return fig


class FocalPlaneGeometryPlot(FocalPlanePlot):
    """Plots the focal plane distribution of a parameter in afw camera
    geometry units: amplifiers and detectors.

    Given the detector positions in x and y, the focal plane positions
    are calculated using the camera model. A 2d binned statistic
    (default is mean) is then calculated and plotted for the parameter
    z as a function of the camera geometry segment the input points
    fall upon.

    The ``xAxisLabel``, ``yAxisLabel``, ``zAxisLabel``, and
    ``statistic`` variables are inherited from the parent class.
    """

    level = ChoiceField[str](
        doc="Which geometry level should values be plotted?",
        default="amplifier",
        allowed={
            "amplifier": "Plot values per readout amplifier.",
            "detector": "Plot values per detector.",
        },
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        base = []
        base.append(("detector", Vector))
        if self.level == "amplifier":
            base.append(("amplifier", Vector))
        base.append(("z", Vector))

        return base

    def makePlot(
        self,
        data: KeyedData,
        camera: Camera,
        plotInfo: Optional[Mapping[str, str]] = None,
        **kwargs,
    ) -> Figure:
        """Prep the catalogue and then make a focalPlanePlot of the given
        column.

        Uses the axisLabels config options `x` and `y` to make an image, where
        the color corresponds to the 2d binned statistic (the mean is the
        default) applied to the `z` column. A summary panel is shown in the
        upper right corner of the resultant plot. The code uses the
        selectorActions to decide which points to plot and the
        statisticSelector actions to determine which points to use for the
        printed statistics.

        Parameters
        ----------
        data : `pandas.core.frame.DataFrame`
            The catalog to plot the points from.  This is expected to
            have the following columns/keys:

            ``"detector"``
                The integer detector id for the points.
            ``"amplifier"``
                The string amplifier name for the points.
            ``"z"``
                The numerical value that will be combined via
                ``statistic`` to the binned value.
            ``"x"``
                Focal plane x position, optional.
            ``"y"``
                Focal plane y position, optional.
        camera : `lsst.afw.cameraGeom.Camera`
            The camera used to map from pixel to focal plane positions.
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:

            ``"run"``
                The output run for the plots (`str`).
            ``"skymap"``
                The type of skymap used for the data (`str`).
            ``"filter"``
                The filter used for this data (`str`).
            ``"tract"``
                The tract that the data comes from (`str`).
            ``"bands"``
                The band(s) that the data comes from (`list` of `str`).

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.
        """

        cmap = mkColormap(["midnightBlue", "lightcyan", "darkgreen"])
        cmap.set_bad(color="none")

        if plotInfo is None:
            plotInfo = {}

        if len(data["z"]) == 0:
            noDataFig = Figure()
            noDataFig.text(0.3, 0.5, "No data to plot after selectors applied")
            noDataFig = addPlotInfo(noDataFig, plotInfo)
            return noDataFig

        if self.addHistogram:
            fig, [ax, histAx] = plt.subplots(1, 2, dpi=300, figsize=(12, 6), width_ratios=[3, 2])
        else:
            fig = plt.figure(dpi=300)
            ax = fig.add_subplot(111)

        detectorIds = np.unique(data["detector"])
        focalPlane_x = np.zeros(len(data["z"]))
        focalPlane_y = np.zeros(len(data["z"]))

        patches = []
        values = []

        # Plot bounding box that will be used to set the axes below.
        plotLimit_x = [0.0, 0.0]
        plotLimit_y = [0.0, 0.0]

        for detectorId in detectorIds:
            detector = camera[detectorId]

            # We can go stright to fp coordinates.
            corners = [(c.getX(), c.getY()) for c in detector.getCorners(FOCAL_PLANE)]
            corners = np.array(corners)

            # U/V coordinates represent focal plane locations.
            minU, minV = corners.min(axis=0)
            maxU, maxV = corners.max(axis=0)

            # See if the plot bounding box needs to be extended:
            if minU < plotLimit_x[0]:
                plotLimit_x[0] = minU
            if minV < plotLimit_y[0]:
                plotLimit_y[0] = minV
            if maxU > plotLimit_x[1]:
                plotLimit_x[1] = maxU
            if maxV > plotLimit_y[1]:
                plotLimit_y[1] = maxV

            # X/Y coordinates represent detector internal coordinates.
            # Detector extent in detector coordinates
            minX, minY = detector.getBBox().getMin()
            maxX, maxY = detector.getBBox().getMax()

            if self.level.lower() == "detector":
                detectorInd = data["detector"] == detectorId

                # This does the appropriate statistic for this
                # detector's data.
                statistic, _, _ = binned_statistic_dd(
                    [focalPlane_x[detectorInd], focalPlane_y[detectorInd]],
                    data["z"][detectorInd],
                    statistic=self.statistic,
                    bins=[1, 1],
                )
                patches.append(Polygon(corners, closed=True))
                values.append(statistic.ravel()[0])
            else:
                # It's at amplifier level.  This uses the focal
                # plane position of the corners of the detector to
                # generate corners for the individual amplifier
                # segments.
                rotation = detector.getOrientation().getNQuarter()  # N * 90 degrees.
                alpha, beta = np.cos(rotation * np.pi / 2.0), np.sin(rotation * np.pi / 2.0)

                # Calculate the rotation matrix between X/Y and U/V
                # coordinates.
                scaleUX = alpha * (maxU - minU) / (maxX - minX)
                scaleVX = beta * (maxV - minV) / (maxX - minX)
                scaleVY = alpha * (maxV - minV) / (maxY - minY)
                scaleUY = beta * (maxU - minU) / (maxY - minY)

                # After the rotation, some of the corners may have
                # negative offsets.  This corresponds to corners that
                # reference the maximum edges of the box in U/V
                # coordinates.
                baseU = minU if rotation % 4 in (0, 1) else maxU
                baseV = maxV if rotation % 4 in (2, 3) else minV

                for amplifier in detector:
                    ampName = amplifier.getName()
                    detectorInd = data["detector"] == detectorId
                    ampInd = data["amplifier"] == ampName
                    ampInd &= detectorInd

                    # Determine amplifier extent in X/Y coordinates.
                    ampMinX, ampMinY = amplifier.getBBox().getMin()
                    ampMaxX, ampMaxY = amplifier.getBBox().getMax()

                    # The corners are rotated into U/V coordinates,
                    # and the appropriate offset added.
                    ampCorners = []
                    ampCorners.append(
                        (
                            scaleUX * (ampMinX - minX) + scaleUY * (ampMinY - minY) + baseU,
                            scaleVY * (ampMinY - minY) + scaleVX * (ampMinX - minX) + baseV,
                        )
                    )
                    ampCorners.append(
                        (
                            scaleUX * (ampMaxX - minX) + scaleUY * (ampMaxY - minY) + baseU,
                            scaleVY * (ampMinY - minY) + scaleVX * (ampMinX - minX) + baseV,
                        )
                    )
                    ampCorners.append(
                        (
                            scaleUX * (ampMaxX - minX) + scaleUY * (ampMaxY - minY) + baseU,
                            scaleVY * (ampMaxY - minY) + scaleVX * (ampMaxX - minX) + baseV,
                        )
                    )
                    ampCorners.append(
                        (
                            scaleUX * (ampMinX - minX) + scaleUY * (ampMinY - minY) + baseU,
                            scaleVY * (ampMaxY - minY) + scaleVX * (ampMaxX - minX) + baseV,
                        )
                    )
                    patches.append(Polygon(ampCorners, closed=True))
                    # This does the appropriate statistic for this
                    # amplifier's data.
                    if len(data["z"][ampInd]) > 0:
                        statistic, _, _ = binned_statistic_dd(
                            [focalPlane_x[ampInd], focalPlane_y[ampInd]],
                            data["z"][ampInd],
                            statistic=self.statistic,
                            bins=[1, 1],
                        )
                        values.append(statistic.ravel()[0])
                    else:
                        values.append(np.nan)

        # Set bounding box for this figure.
        ax.set_xlim(plotLimit_x)
        ax.set_ylim(plotLimit_y)

        # Do not mask values.
        if self.showStats:
            statMed, statMad, statsText = self.statsAndText(values, mask=None)
            bbox = dict(facecolor="paleturquoise", alpha=0.5, edgecolor="none")
            ax.text(0.8, 0.91, statsText, transform=fig.transFigure, fontsize=8, bbox=bbox)

        # Defaults to med + 4 sigma Mad to match
        # the camera team plots
        if self.plotMin is not None:
            vmin = self.plotMin
        else:
            vmin = statMed - 4.0 * statMad
        if self.plotMax is not None:
            vmax = self.plotMax
        else:
            vmax = statMed + 4.0 * statMad

        valuesPlot = np.clip(values, vmin, vmax)

        patchCollection = PatchCollection(
            patches, edgecolor="white", cmap=cmap, linewidth=0.5, linestyle=(0, (0.5, 3))
        )
        patchCollection.set_array(valuesPlot)
        ax.add_collection(patchCollection)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(patchCollection, cax=cax, extend="both")
        text = cax.text(
            0.5,
            0.5,
            self.zAxisLabel,
            color="k",
            rotation="vertical",
            transform=cax.transAxes,
            ha="center",
            va="center",
            fontsize=10,
        )
        text.set_path_effects([pathEffects.Stroke(linewidth=3, foreground="w"), pathEffects.Normal()])
        cax.tick_params(labelsize=7)

        ax.set_xlabel(self.xAxisLabel)
        ax.set_ylabel(self.yAxisLabel)
        ax.tick_params(axis="x", labelrotation=25)
        ax.tick_params(labelsize=7)

        ax.set_aspect("equal")

        if self.addHistogram:
            self._addHistogram(histAx, data["z"])

        plt.draw()

        # Add useful information to the plot
        fig.subplots_adjust(left=0.05, right=0.95)
        fig = plt.gcf()
        if plotInfo:
            fig = addPlotInfo(fig, plotInfo)

        return fig
