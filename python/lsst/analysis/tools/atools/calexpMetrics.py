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

__all__ = (
    "CalexpSummaryMetrics",
    "CalexpMetricsPlot",
)

from ..actions.plot import HistPanel, HistPlot
from ..actions.vector import LoadVector, BandSelector
from ..interfaces import AnalysisTool


class CalexpSummaryMetrics(AnalysisTool):
    """
    Class to load statistics from the summary stats contained with a calexp's
    metadata and write them to metrics.
    """

    propagateData: bool = True

    # raCorners and decCorners statistics cannot be written to a metric,
    # as metrics can only be single-valued (i.e., scalars).
    # Units in comments are to indicate compound units, which are currently
    # unsupported.
    _units = {
        "psfSigma": "pixel",
        "psfArea": "",  # pixel**2
        "psfIxx": "",  # pixel**2
        "psfIyy": "",  # pixel**2
        "psfIxy": "",  # pixel**2
        "ra": "degree",
        "dec": "degree",
        "pixelScale": "",  # arcsec/pixel.
        "zenithDistance": "degree",
        "expTime": "s",
        "zeroPoint": "mag",
        "skyBg": "electron",
        "skyNoise": "electron",
        "meanVar": "",  # electron**2
        "astromOffsetMean": "arcsec",
        "astromOffsetStd": "arcsec",
        "nPsfStar": "ct",
        "psfStarDeltaE1Median": "",
        "psfStarDeltaE2Median": "",
        "psfStarDeltaE1Scatter": "",
        "psfStarDeltaE2Scatter": "",
        "psfStarDeltaSizeMedian": "pixel",
        "psfStarDeltaSizeScatter": "pixel",
        "psfStarScaledDeltaSizeScatter": "",
        "psfTraceRadiusDelta": "pixel",
        "psfApFluxDelta": "",
        "psfApCorrSigmaScaledDelta": "",
        "maxDistToNearestPsf": "pixel",
        "effTime": "s",
        "effTimePsfSigmaScale": "",
        "effTimeSkyBgScale": "",
        "effTimeZeroPointScale": "",
    }

    def setDefaults(self):
        super().setDefaults()

        self.prep.keysToLoad = list(self._units.keys())
        self.produce.metric.units = self._units


class CalexpMetricsPlot(AnalysisTool):

    def setDefaults(self):
        self.prep.selectors.bandSelector = BandSelector()

        self.process.buildActions.posFootprintCount = LoadVector()
        self.process.buildActions.posFootprintCount.vectorKey = (
            "calexpMetadataMetrics_positive_footprint_count"
        )

        self.produce.plot = HistPlot()
        self.produce.plot.panels["panel_footprint_count"] = HistPanel()
        self.produce.plot.panels["panel_footprint_count"].hists = dict(
            posFootprintCount="Number of calexps"
        )
        self.produce.plot.panels["panel_footprint_count"].label = "Positive Footprint Count"
