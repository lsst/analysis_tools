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

from typing import Optional

import lsst.pex.config as pexConfig

from ..actions.vector.calcBinnedStats import CalcBinnedStatsAction
from ..actions.vector.selectors import RangeSelector
from ..analysisParts.diffMatched import (
    MatchedRefCoaddDiffMagTool,
    MatchedRefCoaddDiffPositionTool,
    MatchedRefCoaddTool,
)
from ..interfaces import AnalysisMetric, KeyedData


class MatchedRefCoaddMetric(MatchedRefCoaddTool, AnalysisMetric):
    """Base tool for matched-to-reference metrics on coadds."""

    mag_low_min: int = 15
    mag_low_max: int = 27
    mag_interval: int = 1

    name_prefix = pexConfig.Field[str](default=None, doc="Prefix for metric key")

    names = ("stars", "galaxies", "all")
    types = ("unresolved", "resolved", "all")

    unit = pexConfig.Field[str](default=None, doc="Astropy unit of y-axis values")

    def _validate(self):
        if self.name_prefix is None or self.unit is None:
            raise ValueError(
                f"{self.name_prefix=} and {self.unit=} must not be None;"
                f" did you forget to set a valid context?"
            )

    def configureMetrics(
        self,
        unit: Optional[str] = None,
        name_prefix: Optional[str] = None,
        name_suffix: str = "_mad_ref_mag{minimum}",
        unit_select: str = "mag",
    ):
        """Configure metric actions and return units.

        Parameters
        ----------
        unit : `str`
            The (astropy) unit of the summary statistic metrics.
        name_prefix : `str`
            The prefix for the action (column) name.
        name_suffix : `str`
            The sufffix for the action (column) name.
        unit_select : `str`
            The (astropy) unit of the selection (x-axis) column. Default "mag".

        Returns
        -------
        units : `dict` [`str`, `str`]
            A dict of the unit (value) for each metric name (key)
        """
        unit_is_none = unit is None
        name_prefix_is_none = name_prefix is None

        if unit_is_none or name_prefix_is_none:
            if unit_is_none:
                unit = self.unit
            if name_prefix_is_none:
                name_prefix = self.name_prefix
            self._validate()
        if unit_select is None:
            unit_select = "mag"

        assert name_prefix is not None
        units = {}
        for name, name_class in zip(self.names, self.types):
            name_capital = name.capitalize()
            x_key = f"x{name_capital}"

            for minimum in range(self.mag_low_min, self.mag_low_max + 1):
                action = getattr(self.process.calculateActions, f"{name}{minimum}")
                action.selector_range = RangeSelector(
                    key=x_key,
                    minimum=minimum,
                    maximum=minimum + self.mag_interval,
                )

                action.name_prefix = name_prefix.format(name_class=name_class)
                action.name_suffix = name_suffix.format(minimum=minimum)

                units.update(
                    {
                        action.name_median: unit,
                        action.name_sigmaMad: unit,
                        action.name_count: "count",
                        action.name_select_minimum: unit_select,
                        action.name_select_median: unit_select,
                        action.name_select_maximum: unit_select,
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
                    CalcBinnedStatsAction(key_vector=f"y{name_capital}"),
                )

    def __call__(self, data: KeyedData, **kwargs):
        self._validate()
        return super().__call__(data=data, **kwargs)


# The diamond inheritance on MatchedRefCoaddTool seems ok
class MatchedRefCoaddCModelFluxMetric(MatchedRefCoaddDiffMagTool, MatchedRefCoaddMetric):
    """Metric for diffs between reference and CModel coadd mags."""

    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.unit = "mmag"
        self.name_prefix = "photom_mag_cModelFlux_{name_class}_diff_"
        self.produce.units = self.configureMetrics()

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.unit = ""
        self.name_prefix = "photom_mag_cModelFlux_{name_class}_chi_"
        self.produce.units = self.configureMetrics()

    def setDefaults(self):
        super().setDefaults()


class MatchedRefCoaddPositionMetric(MatchedRefCoaddDiffPositionTool, MatchedRefCoaddMetric):
    """Metric for diffs between reference and base coadd centroids."""

    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.unit = "pix"
        self.name_prefix = f"astrom_{self.variable}_{{name_class}}_diff_"
        self.produce.units = self.configureMetrics()

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.unit = ""
        self.name_prefix = f"astrom_{self.variable}_{{name_class}}_diff_"
        self.produce.units = self.configureMetrics()

    def setDefaults(self):
        super().setDefaults()
