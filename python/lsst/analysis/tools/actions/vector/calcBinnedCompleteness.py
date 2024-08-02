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

__all__ = ("CalcBinnedCompletenessAction",)

from functools import cached_property

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Scalar, Vector, VectorAction
from .selectors import RangeSelector


class CalcBinnedCompletenessAction(KeyedDataAction):
    name_prefix = Field[str](default="", doc="Field name to append statistic names to")
    name_suffix = Field[str](default="", doc="Field name to append to statistic names")
    selector_type_match = ConfigurableActionField[VectorAction](
        doc="Selector to return objects with matching ref and target type classification"
    )
    selector_range_ref = ConfigurableActionField[RangeSelector](doc="Range selector")
    selector_range_target = ConfigurableActionField[RangeSelector](doc="Range selector")

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield from self.selector_range_ref.getInputSchema()
        yield from self.selector_range_target.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return (
            (self.name_mask_ref, Vector),
            (self.name_mask_target, Vector),
            (self.name_count, Scalar),
            (self.name_count, Scalar),
            (self.name_completeness, Scalar),
            (self.name_completeness_bad_match, Scalar),
            (self.name_completeness_good_match, Scalar),
            (self.name_purity, Scalar),
            (self.name_purity_bad_match, Scalar),
            (self.name_purity_good_match, Scalar),
            ("range_maximum", Scalar),
            ("range_minimum", Scalar),
        )

    @cached_property
    def name_count(self):
        return f"{self.name_prefix}count{self.name_suffix}"

    @cached_property
    def name_count_ref(self):
        return f"{self.name_prefix}count_ref{self.name_suffix}"

    @cached_property
    def name_count_target(self):
        return f"{self.name_prefix}count_target{self.name_suffix}"

    @cached_property
    def name_mask_ref(self):
        return f"{self.name_prefix}mask_ref{self.name_suffix}"

    @cached_property
    def name_mask_target(self):
        return f"{self.name_prefix}mask_ref{self.name_suffix}"

    @cached_property
    def name_completeness(self):
        return f"{self.name_prefix}completeness{self.name_suffix}"

    @cached_property
    def name_completeness_bad_match(self):
        return f"{self.name_prefix}completeness_bad_match{self.name_suffix}"

    @cached_property
    def name_completeness_good_match(self):
        return f"{self.name_prefix}completeness_good_match{self.name_suffix}"

    @cached_property
    def name_purity(self):
        return f"{self.name_prefix}purity{self.name_suffix}"

    @cached_property
    def name_purity_bad_match(self):
        return f"{self.name_prefix}purity_bad_match{self.name_suffix}"

    @cached_property
    def name_purity_good_match(self):
        return f"{self.name_prefix}purity_good_match{self.name_suffix}"

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        has_band = "band" in kwargs
        kwargs_format = {}
        if has_band:
            kwargs_format["band"] = kwargs["band"]
        prefix_band = f"{kwargs['band']}_" if has_band else ""

        results = {}
        mask_ref = self.selector_range_ref(data, **kwargs)
        mask_target = self.selector_range_target(data, **kwargs)
        results[self.name_mask_ref.format(**kwargs_format)] = mask_ref
        results[self.name_mask_target.format(**kwargs_format)] = mask_target
        kwargs["mask_ref"] = mask_ref
        kwargs["mask_target"] = mask_target

        n_ref = np.sum(mask_ref)
        n_target = np.sum(mask_target)
        mask_any = mask_ref | mask_target
        matched = mask_ref & mask_target

        n_matched = np.sum(matched)
        n_matched_same = np.sum(data[self.key_class_ref][matched] == data[self.key_class_target][matched])
        n_matched_diff = n_matched - n_matched_same

        results[self.name_count.format(**kwargs_format)] = np.sum(mask_any)
        results[self.name_count_ref.format(**kwargs_format)] = n_ref
        results[self.name_count_target.format(**kwargs_format)] = n_target
        results[self.name_completeness.format(**kwargs_format)] = n_matched / n_ref
        results[self.name_completeness_bad_match.format(**kwargs_format)] = n_matched_diff / n_ref
        results[self.name_completeness_good_match.format(**kwargs_format)] = n_matched_same / n_ref
        results[self.name_purity.format(**kwargs_format)] = n_matched / n_target
        results[self.name_purity_bad_match.format(**kwargs_format)] = n_matched_diff / n_target
        results[self.name_purity_good_match.format(**kwargs_format)] = n_matched_same / n_target

        results[f"{prefix_band}range_maximum"] = self.selector_range.maximum
        results[f"{prefix_band}range_minimum"] = self.selector_range.minimum

        return results
