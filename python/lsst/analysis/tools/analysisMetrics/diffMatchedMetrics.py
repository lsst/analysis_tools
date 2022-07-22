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
    mag_low_min: int = 15
    mag_low_max: int = 27
    mag_interval: int = 1

    names = ("stars", "galaxies")
    types = ("unresolved", "resolved")

    def setDefaults(self):
        super().setDefaults()

        for name in self.names:
            name_capital = name.capitalize()
            for minimum in range(self.mag_low_min, self.mag_low_max + 1):
                setattr(
                    self.process.calculateActions,
                    f"{name}{minimum}",
                    CalcBinnedStatsAction(vectorKey=f"y{name_capital}"),
                )


class MatchedRefCoaddDiffCModelFluxMetric(MatchedRefCoaddDiffMagTool, MatchedRefCoaddMetric):
    def matchedRefDiffMagContext(self):
        super(MatchedRefCoaddDiffCModelFluxMetric, self).matchedRefDiffMagContext()

    def matchedRefDiffFluxChiContext(self):
        super(MatchedRefCoaddDiffCModelFluxMetric, self).matchedRefDiffFluxChiContext()

    def setDefaults(self):
        super(MatchedRefCoaddDiffCModelFluxMetric, self).setDefaults()

        units = {}

        for name, name_class in zip(self.names, self.types):
            name_capital = name.capitalize()
            x_key = f"x{name_capital}"

            for minimum in range(self.mag_low_min, self.mag_low_max + 1):
                action = getattr(self.process.calculateActions, f"{name}{minimum}")
                name_prefix = f"photom_mag_cModelFlux_{name_class}_diff_sig_mad_ref_mag{minimum}"
                action.rangeSelector = RangeSelector(
                    column=x_key,
                    minimum=minimum,
                    maximum=minimum + self.mag_interval,
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
