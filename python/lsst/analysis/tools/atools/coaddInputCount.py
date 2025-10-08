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

__all__ = ("CoaddInputCount", "CoaddQualityCheck", "CoaddQualityPlot")

from lsst.pex.config import ListField

from ..actions.plot.calculateRange import MinMax
from ..actions.plot.coaddDepthPlot import CoaddDepthPlot
from ..actions.plot.skyPlot import SkyPlot
from ..actions.scalar.scalarActions import MeanAction, MedianAction, SigmaMadAction, StdevAction
from ..actions.vector import BandSelector, CoaddPlotFlagSelector, DownselectVector, LoadVector, SnSelector
from ..interfaces import AnalysisTool


class CoaddInputCount(AnalysisTool):
    """Tract-wide metrics pertaining to how many exposures have gone into
    a deep coadd, per band.

    This AnalysisTool is designed to run on an object table, which is only
    created for deep coadds, not template coadds.
    """

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        # Set this to an empty list to look at the band
        # the plot is being made in.
        self.prep.selectors.flagSelector.bands = []

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 5

        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "coord_ra"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "coord_dec"
        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.fluxType = "{band}_psfFlux"
        self.process.buildActions.statMask.threshold = 5

        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "{band}_inputCount"

        self.process.buildActions.patch = LoadVector()
        self.process.buildActions.patch.vectorKey = "patch"

        self.process.calculateActions.median = MedianAction()
        self.process.calculateActions.median.vectorKey = "z"

        self.process.calculateActions.mean = MeanAction()
        self.process.calculateActions.mean.vectorKey = "z"

        self.process.calculateActions.sigmaMad = SigmaMadAction()
        self.process.calculateActions.sigmaMad.vectorKey = "z"

        # SkyPlot of number of contributing exposures in coad, per tract/band:
        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["any"]
        self.produce.plot.plotName = "{band}_inputCount"
        self.produce.plot.xAxisLabel = "R.A. (deg)"
        self.produce.plot.yAxisLabel = "Dec. (deg)"
        self.produce.plot.zAxisLabel = "Input Count"
        self.produce.plot.plotOutlines = True
        self.produce.plot.showExtremeOutliers = False
        self.produce.plot.colorbarRange = MinMax()

        # Summary metrics for the whole coadd, per tract/band.
        self.produce.metric.units = {"median": "ct", "sigmaMad": "ct", "mean": "ct"}
        self.produce.metric.newNames = {
            "median": "{band}_inputCount_median",
            "mean": "{band}_inputCount_mean",
            "sigmaMad": "{band}_inputCount_sigmaMad",
        }


class CoaddQualityCheck(AnalysisTool):
    """Compute the percentage of each coadd that has a number of input
    exposures exceeding a threshold.

    This AnalysisTool is designed to run on any coadd, provided a
    coadd_depth_table is created first (via CoaddDepthSummaryAnalysisTask).

    For example, if exactly half of a coadd patch contains 15 overlapping
    constituent visits and half contains fewer, the value computed for
    `depth_above_threshold_12` would be 50.

    These values come from the n_image data product, which is an image
    identical to the coadd but with pixel values of the number of input
    images instead of flux or counts. n_images are persisted during
    coadd assembly.
    """

    band_list = ListField(
        default=["u", "g", "r", "i", "z", "y"],
        dtype=str,
        doc="Bands for colors.",
    )

    threshold_list = ListField(
        default=[1, 3, 5, 12],
        dtype=int,
        doc="The n_image pixel value thresholds.",
    )

    quantile_list = ListField(
        default=[5, 10, 25, 50, 75, 90, 95],
        dtype=int,
        doc="The percentiles at which to compute n_image values, in ascending order.",
    )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.patch = LoadVector(vectorKey="patch")
        self.process.buildActions.band = LoadVector(vectorKey="band")

    def finalize(self):
        for threshold in self.threshold_list:
            # This gives a RuntimeWarning whenever the tract doesn't have
            # a particular band. Need to use a band list derived from the
            # bands found in the "band" column of a given tract.
            for band in self.band_list:
                name = f"depth_above_threshold_{threshold}"
                setattr(self.process.buildActions, name, LoadVector(vectorKey=name))
                setattr(
                    self.process.filterActions,
                    f"{name}_{band}",
                    DownselectVector(vectorKey=name, selector=BandSelector(bands=[band])),
                )
                setattr(
                    self.process.calculateActions,
                    f"{name}_{band}_median",
                    MedianAction(vectorKey=f"{name}_{band}"),
                )
                setattr(
                    self.process.calculateActions,
                    f"{name}_{band}_mean",
                    MeanAction(vectorKey=f"{name}_{band}"),
                )
                setattr(
                    self.process.calculateActions,
                    f"{name}_{band}_stdev",
                    StdevAction(vectorKey=f"{name}_{band}"),
                )

                # The units for the quantity are dimensionless (percentage)
                self.produce.metric.units[f"{name}_{band}_median"] = ""
                self.produce.metric.units[f"{name}_{band}_mean"] = ""
                self.produce.metric.units[f"{name}_{band}_stdev"] = ""

        for quantile in self.quantile_list:
            for band in self.band_list:
                name = f"depth_{quantile}_percentile"
                setattr(self.process.buildActions, name, LoadVector(vectorKey=name))
                setattr(
                    self.process.filterActions,
                    f"{name}_{band}",
                    DownselectVector(vectorKey=name, selector=BandSelector(bands=[band])),
                )
                setattr(
                    self.process.calculateActions,
                    f"{name}_{band}_median",
                    MedianAction(vectorKey=f"{name}_{band}"),
                )
                setattr(
                    self.process.calculateActions,
                    f"{name}_{band}_mean",
                    MeanAction(vectorKey=f"{name}_{band}"),
                )
                setattr(
                    self.process.calculateActions,
                    f"{name}_{band}_stdev",
                    StdevAction(vectorKey=f"{name}_{band}"),
                )

                # The units for the quantity are dimensionless (percentage)
                self.produce.metric.units[f"{name}_{band}_median"] = ""
                self.produce.metric.units[f"{name}_{band}_mean"] = ""
                self.produce.metric.units[f"{name}_{band}_stdev"] = ""


class CoaddQualityPlot(AnalysisTool):
    """Make a plot of coadd depth."""

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.patch = LoadVector(vectorKey="patch")
        self.process.buildActions.band = LoadVector(vectorKey="band")
        self.process.buildActions.depth = LoadVector(vectorKey="depth")
        self.process.buildActions.pixels = LoadVector(vectorKey="pixels")

        self.produce.plot = CoaddDepthPlot()
