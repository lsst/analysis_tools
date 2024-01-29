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


from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, ScalarType, Vector
from astropy.table import Table, vstack
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
from .plotUtils import addPlotInfo
from typing import Mapping

__all__ = ("PercentilePlot",)


class PercentilePlot(PlotAction):
    """Makes a scatter plot of the data with a marginal
    histogram for each axis.
    """

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | ScalarType]] = []
        base.append(("amplifier", Vector))
        base.append(("detector", Vector))
        base.append(("percentile_0", Vector))
        base.append(("percentile_5", Vector))
        base.append(("percentile_16", Vector))
        base.append(("percentile_50", Vector))
        base.append(("percentile_84", Vector))
        base.append(("percentile_95", Vector))
        base.append(("percentile_100", Vector))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """NOTE currently can only check that something is not a Scalar, not
        check that the data is consistent with Vector
        """
        needed = self.getFormattedInputSchema(**kwargs)
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in input")
        for name, typ in needed:
            isScalar = issubclass((colType := type(data[name.format(**kwargs)])), Scalar)
            if isScalar and typ != Scalar:
                raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}")

    def makePlot(self, data, plotInfo, **kwargs):
        """Makes a plot showing the percentiles of the normalized distribution
        of the data.

        Parameters
        ----------
        data : `KeyedData`
            All the data
        plotInfo : `dict`
            A dictionary of information about the data being plotted with keys:
            ``camera``
                The camera used to take the data (`lsst.afw.cameraGeom.Camera`)
            ``"cameraName"``
                The name of camera used to take the data (`str`).
            ``"filter"``
                The filter used for this data (`str`).
            ``"ccdKey"``
                The ccd/dectector key associated with this camera (`str`).
            ``"visit"``
                The visit of the data; only included if the data is from a
                single epoch dataset (`str`).
            ``"patch"``
                The patch that the data is from; only included if the data is
                from a coadd dataset (`str`).
            ``"tract"``
                The tract that the data comes from (`str`).
            ``"photoCalibDataset"``
                The dataset used for the calibration, e.g. "jointcal" or "fgcm"
                (`str`).
            ``"skyWcsDataset"``
                The sky Wcs dataset used (`str`).
            ``"rerun"``
                The rerun the data is stored in (`str`).

        Returns
        ------
        ``fig``
            The figure to be saved (`matplotlib.figure.Figure`).

        Notes
        -----
        Makes a plot showing the normalized percentile distribution of data.
        """
        amplifiers = [
            "C17",
            "C07",
            "C16",
            "C06",
            "C15",
            "C05",
            "C14",
            "C04",
            "C13",
            "C03",
            "C12",
            "C02",
            "C11",
            "C01",
            "C10",
            "C00",
        ]
        # TODO: generalize to make N per-detector plots
        detector = data["detector"] == 0
        data = vstack([Table(data)[detector & (data["amplifier"] == amp)][0] for amp in amplifiers])
        percentiles = ["0", "5", "16", "50", "84", "95", "100"]
        distributions = [data[f"percentile_{pct}"] for pct in percentiles]
        medians = [np.nanmedian(dist) for dist in distributions]
        normalizedDistributions = [np.abs(dist / med) for (med, dist) in list(zip(medians, distributions))]

        fig, axs = plt.subplots(nrows=8, ncols=2, sharex=True, sharey=True)
        # Set threshold for a hot column.
        threshold = [0.1, 10]
        pcts = np.array([int(pct) for pct in percentiles])
        for i, ax in enumerate(axs.reshape(16)):
            # Get the distribution for a single amplifier.
            distribution = np.array([dist[i] for dist in normalizedDistributions])

            # Plot points below, above, and within the threshold distinctly.
            belowThreshold = np.where(distribution < threshold[0])[0]
            aboveThreshold = np.where(distribution > threshold[1])[0]
            withinThreshold = np.where((distribution > threshold[0]) & (distribution < threshold[1]))
            ax.scatter(
                pcts[belowThreshold],
                distribution[belowThreshold],
                c="r",
                marker="v",
                label="outside threshold" if i == 0 else "",
            )
            ax.scatter(pcts[aboveThreshold], distribution[aboveThreshold], c="r", marker="^")
            ax.scatter(
                pcts[withinThreshold],
                distribution[withinThreshold],
                c="C0",
                marker="o",
                s=10,
                label="within threshold" if i == 0 else "",
            )
            # Connect the scattered dots.
            ax.plot(pcts, distribution, zorder=0)
            # Plot the ideal line.
            ax.hlines(
                1.0, xmin=pcts[0], xmax=pcts[-1], colors="k", linestyle="--", label="1" if i == 0 else ""
            )
            ax.set_ylabel(data["amplifier"][i])
            ax.set_yscale("log")
            ax.tick_params("x", labelrotation=45)

        plt.xticks(ticks=pcts, labels=percentiles)
        fig.supxlabel("Percentile")
        fig.supylabel("Normalized distribution")
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0)
        plt.figlegend()

        # Add useful information to the plot
        fig = plt.gcf()
        addPlotInfo(fig, plotInfo)
        return fig
