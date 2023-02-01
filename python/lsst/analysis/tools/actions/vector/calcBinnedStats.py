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

__all__ = ("CalcBinnedStatsAction",)

from functools import cached_property
from typing import cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Scalar, Vector
from ..keyedData.summaryStatistics import SummaryStatisticAction
from .selectors import RangeSelector


class CalcBinnedStatsAction(KeyedDataAction):
    key_vector = Field[str](doc="Vector on which to compute statistics")
    name_prefix = Field[str](doc="Field name to append stat names to")
    name_suffix = Field[str](doc="Field name to append to stat names")
    selector_range = ConfigurableActionField[RangeSelector](doc="Range selector")

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield (self.key_vector, Vector)
        yield from self.selector_range.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return (
            (self.name_mask, Vector),
            (self.name_median, Scalar),
            (self.name_sig_mad, Scalar),
            (self.name_count, Scalar),
            (self.name_select_maximum, Scalar),
            (self.name_select_median, Scalar),
            (self.name_select_minimum, Scalar),
            ("range_maximum", Scalar),
            ("range_minimum", Scalar),
        )

    @cached_property
    def name_count(self):
        return f"{self.name_prefix}count{self.name_suffix}"

    @cached_property
    def name_mask(self):
        return f"{self.name_prefix}mask{self.name_suffix}"

    @cached_property
    def name_median(self):
        return f"{self.name_prefix}median{self.name_suffix}"

    @cached_property
    def name_select_maximum(self):
        return f"{self.name_prefix}select_maximum{self.name_suffix}"

    @cached_property
    def name_select_median(self):
        return f"{self.name_prefix}select_median{self.name_suffix}"

    @cached_property
    def name_select_minimum(self):
        return f"{self.name_prefix}select_minimum{self.name_suffix}"

    @cached_property
    def name_sigmaMad(self):
        return f"{self.name_prefix}sigmaMad{self.name_suffix}"

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        mask = self.selector_range(data, **kwargs)
        results[self.name_mask] = mask
        kwargs["mask"] = mask

        action = SummaryStatisticAction(vectorKey=self.key_vector)
        # this is sad, but pex_config seems to have broken behavior that
        # is dangerous to fix
        action.setDefaults()

        for name, value in action(data, **kwargs).items():
            results[getattr(self, f"name_{name}")] = value

        values = cast(Vector, data[self.selector_range.key][mask])  # type: ignore
        valid = np.sum(np.isfinite(values)) > 0
        results[self.name_select_maximum] = cast(Scalar, float(np.nanmax(values)) if valid else np.nan)
        results[self.name_select_median] = cast(Scalar, float(np.nanmedian(values)) if valid else np.nan)
        results[self.name_select_minimum] = cast(Scalar, float(np.nanmin(values)) if valid else np.nan)
        results["range_maximum"] = self.selector_range.maximum
        results["range_minimum"] = self.selector_range.minimum

        return results
