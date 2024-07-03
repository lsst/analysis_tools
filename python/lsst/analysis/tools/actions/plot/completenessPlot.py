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


from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
from lsst.pex.config import Field, ListField
from matplotlib.figure import Figure

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, ScalarType, Vector
from .plotUtils import addPlotInfo

__all__ = ("CompletenessHist",)


class CompletenessHist(PlotAction):
    """Makes a scatter plot of the data with a marginal
    histogram for each axis.
    """

    magKey = Field[str](doc="Name of the magnitude column.", default="mag")
    matchDistanceKey = Field[str](doc="Name of the match distance column.", default="matchDistance")
    xAxisLabel = Field[str](doc="Label for the x axis.", default="Input Magnitude (mag)")
    inputLabel = Field[str](doc="Label for the input source histogram.", default="Synthetic Inputs")
    outputLabel = Field[str](doc="Label for the recovered source histogram.", default="Synthetic Recovered")
    numBins = Field[int](doc="Number of bins to use for the histograms.", default=100)
    completenessPercentiles = ListField[float](
        doc="Record the magnitudes at these percentiles", default=[84.0, 50.0, 16.0]
    )

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | ScalarType]] = []
        base.append((self.magKey, Vector))
        base.append((self.matchDistanceKey, Vector))
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
        """Makes a plot showing the fraction of injected sources recovered by
        input magnitude.

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
        Makes a histogram showing the fraction recovered in each magnitude
        bin with the number input and recovered overplotted.

        Examples
        --------
        An example of the plot produced from this code is here:

        .. image:: /_static/analysis_tools/completenessPlotExample.png

        For a detailed example of how to make a plot from the command line
        please see the
        :ref:`getting started guide<analysis-tools-getting-started>`.
        """

        # Make plot showing the fraction recovered in magnitude bins
        fig, axLeft = plt.subplots(dpi=300)
        axLeft.tick_params(axis="y", labelcolor="C0")
        axLeft.set_xlabel(self.xAxisLabel)
        axLeft.set_ylabel("Fraction Recovered", color="C0")
        axRight = axLeft.twinx()
        axRight.set_ylabel("Number of Sources")
        matched = np.isfinite(data[self.matchDistanceKey])
        nInput, bins, _ = axRight.hist(
            data[self.magKey],
            range=(np.nanmin(data[self.magKey]), np.nanmax(data[self.magKey])),
            bins=self.numBins,
            log=True,
            histtype="step",
            label=self.inputLabel,
            color="black",
        )
        nOutput, _, _ = axRight.hist(
            data[self.magKey][matched],
            range=(np.nanmin(data[self.magKey][matched]), np.nanmax(data[self.magKey][matched])),
            bins=bins,
            log=True,
            histtype="step",
            label=self.outputLabel,
            color="grey",
        )

        # Find bin where the fraction recovered falls below a given percentile.
        percentileInfo = []
        xlims = plt.gca().get_xlim()
        for pct in self.completenessPercentiles:
            pct /= 100
            magArray = np.where((nOutput / nInput < pct))[0]
            if len(magArray) == 0:
                mag = np.nan
            else:
                mag = np.min(bins[magArray])
                axLeft.plot([xlims[0], mag], [pct, pct], ls=":", color="grey")
                axLeft.plot([mag, mag], [0, pct], ls=":", color="grey")
                percentileInfo.append("Magnitude at {}% recovered: {:0.2f}".format(pct * 100, mag))
        plt.xlim(xlims)
        axLeft.set_ylim(0, 1.05)
        axLeft.axhline(1, color="grey", ls="--")
        axLeft.bar(
            bins[:-1],
            nOutput / nInput,
            width=np.diff(bins),
            align="edge",
            color="C0",
            alpha=0.5,
            zorder=10,
        )

        # Add useful information to the plot
        fig = plt.gcf()
        addPlotInfo(fig, plotInfo)
        fig.legend()
        statsText = ""
        for info in percentileInfo:
            statsText += f"{info}\n"
        bbox = dict(edgecolor="grey", linestyle=":", facecolor="none")
        fig.text(0.7, 0.075, statsText[:-1], bbox=bbox, transform=fig.transFigure, fontsize=6)
        fig.subplots_adjust(bottom=0.2)
        return fig
