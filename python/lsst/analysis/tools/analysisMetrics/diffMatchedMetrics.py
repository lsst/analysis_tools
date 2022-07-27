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

__all__ = ("MatchedRefCoaddMetric", "MatchedRefCoaddCModelFluxMetric", "MatchedRefCoaddPositionMetric")

from ..actions.vector.calcBinnedStats import CalcBinnedStatsAction
from ..actions.vector.selectors import RangeSelector
from ..analysisParts.diffMatched import MatchedRefCoaddDiffMagTool, MatchedRefCoaddDiffPositionTool
from ..interfaces import AnalysisMetric


class MatchedRefCoaddMetric(AnalysisMetric):
    mag_low_min: int = 15
    mag_low_max: int = 27
    mag_interval: int = 1

    names = ("stars", "galaxies", "all")
    types = ("unresolved", "resolved", "all")

    def configureMetrics(self, unit: str, name_prefix: str, name_suffix: str):
        units = {}
        for name, name_class in zip(self.names, self.types):
            name_capital = name.capitalize()
            x_key = f"x{name_capital}"

            for minimum in range(self.mag_low_min, self.mag_low_max + 1):
                action = getattr(self.process.calculateActions, f"{name}{minimum}")
                action.rangeSelector = RangeSelector(
                    column=x_key,
                    minimum=minimum,
                    maximum=minimum + self.mag_interval,
                )
                action.name_prefix = name_prefix.format(name_class=name_class)
                action.name_suffix = name_suffix.format(minimum=minimum)

                units.update(
                    {
                        action.name_selectMedian: unit,
                        action.name_median: unit,
                        action.name_sigmaMad: unit,
                        action.name_count: "count",
                    }
                )
        return units

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


class MatchedRefCoaddCModelFluxMetric(MatchedRefCoaddDiffMagTool, MatchedRefCoaddMetric):
    def matchedRefDiffContext(self):
        super(MatchedRefCoaddCModelFluxMetric, self).matchedRefDiffContext()
        self.produce.units = self.configureMetrics(
            unit="mag",
            name_prefix="photom_mag_cModelFlux_{name_class}_diff_",
            name_suffix="_mad_ref_mag{minimum}",
        )

    def matchedRefChiContext(self):
        super(MatchedRefCoaddCModelFluxMetric, self).matchedRefChiContext()
        self.produce.units = self.configureMetrics(
            unit="",
            name_prefix="photom_mag_cModelFlux_{name_class}_chi_",
            name_suffix="_mad_ref_mag{minimum}",
        )

    def setDefaults(self):
        super(MatchedRefCoaddCModelFluxMetric, self).setDefaults()


class MatchedRefCoaddPositionMetric(MatchedRefCoaddDiffPositionTool, MatchedRefCoaddMetric):
    def matchedRefDiffContext(self):
        super(MatchedRefCoaddPositionMetric, self).matchedRefDiffContext()
        self.produce.units = self.configureMetrics(
            unit="mag",
            name_prefix=f"astrom_{self.variable}_{{name_class}}_diff_",
            name_suffix="_mad_ref_mag{minimum}",
        )

    def matchedRefChiContext(self):
        super(MatchedRefCoaddPositionMetric, self).matchedRefChiContext()
        self.produce.units = self.configureMetrics(
            unit="",
            name_prefix=f"astrom_{self.variable}_{{name_class}}_diff_",
            name_suffix="_mad_ref_mag{minimum}",
        )

    def setDefaults(self):
        super(MatchedRefCoaddPositionMetric, self).setDefaults()
