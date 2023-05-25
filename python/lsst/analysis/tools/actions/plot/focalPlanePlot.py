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

__all__ = ("FocalPlanePlot",)

from typing import Mapping, Optional

import matplotlib.patheffects as pathEffects
import matplotlib.pyplot as plt
import numpy as np
from lsst.afw.cameraGeom import FOCAL_PLANE, PIXELS, Camera
from lsst.pex.config import Field
from matplotlib.figure import Figure
from scipy.stats import binned_statistic_2d

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from ...statistics import nansigmaMad
from .plotUtils import addPlotInfo, sortAllArrays


class FocalPlanePlot(PlotAction):
    """Plots the focal plane distribution of a parameter.

    Given the detector positions in x and y, the focal plane positions are
    calculated using the camera model. A 2d binned statistic (default is mean)
    is then calculated and plotted for the parameter z as a function of the
    focal plane coordinates.
    """

    xAxisLabel = Field[str](doc="Label to use for the x axis.", optional=False)
    yAxisLabel = Field[str](doc="Label to use for the y axis.", optional=False)
    zAxisLabel = Field[str](doc="Label to use for the z axis.", optional=False)

    nBins = Field[int](
        doc="Number of bins to use within the effective plot ranges along the spatial directions.",
        default=200,
    )
    statistic = Field[str](
        doc="Operation to perform in binned_statistic_2d",
        default="mean",
    )

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
        med = np.nanmedian(arr)
        sigMad = nansigmaMad(arr)

        statsText = (
            "Median: {:0.2f}\n".format(med)
            + r"$\sigma_{MAD}$: "
            + "{:0.2f}\n".format(sigMad)
            + r"n$_{points}$: "
            + "{}".format(numPoints)
        )

        return med, sigMad, statsText

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

        if len(data["x"]) == 0:
            noDataFig = Figure()
            noDataFig.text(0.3, 0.5, "No data to plot after selectors applied")
            noDataFig = addPlotInfo(noDataFig, plotInfo)
            return noDataFig

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

        # Add an arbitrary small offset to bins to ensure that the minimum does
        # not equal the maximum.
        binsx = np.linspace(focalPlane_x.min() - 1e-5, focalPlane_x.max() + 1e-5, self.nBins)
        binsy = np.linspace(focalPlane_y.min() - 1e-5, focalPlane_y.max() + 1e-5, self.nBins)

        statistic, x_edge, y_edge, binnumber = binned_statistic_2d(
            focalPlane_x, focalPlane_y, data["z"], statistic=self.statistic, bins=[binsx, binsy]
        )
        binExtent = [x_edge[0], x_edge[-1], y_edge[0], y_edge[-1]]

        sortedArrs = sortAllArrays([data["z"], data["x"], data["y"], data["statMask"]])
        [colorVals, xs, ys, stat] = sortedArrs
        statMed, statMad, statsText = self.statsAndText(colorVals, mask=stat)
        bbox = dict(facecolor="paleturquoise", alpha=0.5, edgecolor="none")
        ax.text(0.8, 0.91, statsText, transform=fig.transFigure, fontsize=8, bbox=bbox)

        median = np.nanmedian(statistic.ravel())
        mad = nansigmaMad(statistic.ravel())

        vmin = median - 2 * mad
        vmax = median + 2 * mad

        plot = ax.imshow(statistic.T, extent=binExtent, vmin=vmin, vmax=vmax, origin="lower")

        cax = fig.add_axes([0.87 + 0.04, 0.11, 0.04, 0.77])
        plt.colorbar(plot, cax=cax, extend="both")
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
        plt.draw()

        # Add useful information to the plot
        plt.subplots_adjust(wspace=0.0, hspace=0.0, right=0.85)
        fig = plt.gcf()
        fig = addPlotInfo(fig, plotInfo)

        return fig
