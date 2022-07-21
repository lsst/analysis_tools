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

from ..actions.vector.calcBinnedStats import CalcBinnedStatsAction
from ..actions.vector.selectors import RangeSelector
from ..analysisParts.diffMatched import MatchedRefCoaddDiffMagTool
from ..interfaces import AnalysisMetric


class MatchedRefCoaddMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()

        self.process.calculateActions.galaxies = CalcBinnedStatsAction(vectorKey="yGalaxies")
        self.process.calculateActions.stars = CalcBinnedStatsAction(vectorKey="yStars")


class MatchedRefCoaddDiffCModelFluxMetric(MatchedRefCoaddDiffMagTool, MatchedRefCoaddMetric):
    def setDefaults(self):
        super(MatchedRefCoaddDiffCModelFluxMetric).setDefaults()
        minimum = 15
        maximum = 16

        units = {}

        for action, x_key, name_class in (
            (self.process.calculateActions.galaxies, "xGalaxies", "resolved"),
            (self.process.calculateActions.stars, "xStars", "unresolved"),
        ):
            name_prefix = f"photom_mag_cModelFlux_{name_class}_diff_sig_mad_ref_mag15"
            action.rangeSelector = RangeSelector(
                column=x_key,
                minimum=minimum,
                maximum=maximum,
            )
            action.name_prefix = name_prefix

            units.update(
                {
                    action.name_selectMedian: "mag",
                    action.name_median: "mag",
                    action.name_sigmaMad: "mag",
                    action.name_count: "count",
                }
            )

        self.produce.units = units  # type: ignore
