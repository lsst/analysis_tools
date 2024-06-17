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

__all__ = ("CalexpSummaryMetrics",)

from ..interfaces import AnalysisTool


class CalexpSummaryMetrics(AnalysisTool):
    """
    Class to load statistics from the summary stats contained with a calexp's
    metadata and write them to metrics.
    """

    propagateData: bool = True

    # raCorners and decCorners statistics cannot be written to a metric,
    # as metrics can only be single-valued (i.e., scalars).
    _units = {
        "psfSigma": "pixel",
        "psfArea": "",  # There is no astropy sq. pixel unit.
        "psfIxx": "pixel",
        "psfIyy": "pixel",
        "psfIxy": "pixel",
        "ra": "degree",
        "dec": "degree",
        "zenithDistance": "degree",
        "expTime": "s",
        "zeroPoint": "mag",
        "skyBg": "mag",
        "skyNoise": "mag",
        "meanVar": "mag",
        "astromOffsetMean": "arcsec",
        "astromOffsetStd": "arcsec",
        "nPsfStar": "",
        "psfStarDeltaE1Median": "pixel",
        "psfStarDeltaE2Median": "pixel",
        "psfStarDeltaE1Scatter": "pixel",
        "psfStarDeltaE2Scatter": "pixel",
        "psfStarDeltaSizeMedian": "pixel",
        "psfStarDeltaSizeScatter": "pixel",
        "psfStarScaledDeltaSizeScatter": "pixel",
        "psfTraceRadiusDelta": "pixel",
        "psfApFluxDelta": "",  # This is a normalized-to-one value.
        "psfApCorrSigmaScaledDelta": "",  # This is a multiplicative factor.
        "maxDistToNearestPsf": "deg",
        "effTime": "s",
        "effTimePsfSigmaScale": "s",
        "effTimeSkyBgScale": "s",
        "effTimeZeroPointScale": "s",
    }

    def setDefaults(self):
        super().setDefaults()

        self.prep.keysToLoad = list(self._units.keys())
        self.produce.metric.units = self._units
