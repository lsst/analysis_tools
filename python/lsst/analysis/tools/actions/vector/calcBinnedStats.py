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

import numpy as np
from lsst.pex.config import Field
from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Scalar, Vector
from ..keyedData.summaryStatistics import SummaryStatisticAction
from .selectors import RangeSelector


class CalcBinnedStatsAction(KeyedDataAction):
    name_prefix = Field[str](doc="Field name to append stat names to", default="")
    name_suffix = Field[str](doc="Field name to append to stat names", default="")
    rangeSelector = ConfigurableActionField[RangeSelector](doc="Range selector")
    vectorKey = Field[str](doc="Vector on which to compute statistics")

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield (self.vectorKey, Vector)
        yield from self.rangeSelector.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return (
            (self.name_mask, Vector),
            (self.name_median, Scalar),
            (self.name_sig_mad, Scalar),
            (self.name_count, Scalar),
            (self.name_select_median, Scalar),
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
    def name_selectMedian(self):
        return f"{self.name_prefix}selectMedian{self.name_suffix}"

    @cached_property
    def name_sigmaMad(self):
        return f"{self.name_prefix}sigmaMad{self.name_suffix}"

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        mask = self.rangeSelector(data, **kwargs)
        results[self.name_mask] = mask
        kwargs["mask"] = mask

        statAction = SummaryStatisticAction(vectorKey=self.vectorKey)
        # this is sad, but pex_config seems to have broken behavior that
        # is dangerous to fix
        statAction.setDefaults()

        for name, value in statAction(data, **kwargs).items():
            results[getattr(self, f"name_{name}")] = value

        results[self.name_selectMedian] = float(np.nanmedian(data[self.rangeSelector.column][mask]))
        results["range_maximum"] = self.rangeSelector.maximum  # type: ignore
        results["range_minimum"] = self.rangeSelector.minimum  # type: ignore

        return results
