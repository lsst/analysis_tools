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

__all__ = ("MatchedRefCoaddMetric",)

from ..analysisParts.diffMatched import setMatchedRefCoaddDefaults, setMatchedRefCoaddDiffMagDefaults
from ..interfaces import AnalysisMetric
from ..actions.vector.calcBinnedStats import CalcBinnedStatsAction
from ..actions.vector.selectors import RangeSelector


class MatchedRefCoaddMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()
        setMatchedRefCoaddDefaults(self)

        self.process.calculateActions.galaxies = CalcBinnedStatsAction(vectorKey="yGalaxies")
        self.process.calculateActions.stars = CalcBinnedStatsAction(vectorKey="yStars")

        # TODO: Finish
        self.produce.units = {  # type: ignore
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        }


class MatchedRefCoaddDiffCModelFluxMetric(MatchedRefCoaddMetric):
    def setDefaults(self):
        super().setDefaults()
        setMatchedRefCoaddDiffMagDefaults(self)
        minimum = 15
        maximum = 16

        for action, name_class in (
            (self.process.calculateActions.galaxies, "resolved"),
            (self.process.calculateActions.stars, "unresolved"),
        ):
            name_prefix = f"photom_mag_cModelFlux_{name_class}_diff_sig_mad_ref_mag15"
            action.rangeSelector = RangeSelector(
                column="mags_ref",
                minimum=minimum,
                maximum=maximum,
            )
            action.name_prefix = name_prefix
