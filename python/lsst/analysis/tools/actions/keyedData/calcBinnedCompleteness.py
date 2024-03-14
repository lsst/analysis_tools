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

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Scalar, Vector
from ...math import divide
from ..vector.selectors import RangeSelector


class CalcBinnedCompletenessAction(KeyedDataAction):
    key_match_distance = Field[str](
        default="match_distance",
        doc="Key for column with distance between matched objects",
    )
    key_matched_class = Field[str](
        default="matched_class",
        doc="Key for boolean vector (True if matched objects have the same class as their ref match)",
    )
    key_mask_ref = Field[str](
        default=None,
        doc="Key for mask to apply for reference objects in completeness",
        optional=True,
    )
    key_mask_target = Field[str](
        default=None,
        doc="Key for mask to apply for target objects in purity",
        optional=True,
    )
    name_prefix = Field[str](default="", doc="Field name to append statistic names to")
    name_suffix = Field[str](default="", doc="Field name to append to statistic names")
    selector_range_ref = ConfigurableActionField[RangeSelector](doc="Range selector")
    selector_range_target = ConfigurableActionField[RangeSelector](doc="Range selector")

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield self.key_match_distance, Vector
        yield self.key_matched_class, Vector
        if self.key_mask_ref:
            yield self.key_mask_ref, Vector
        if self.key_mask_target:
            yield self.key_mask_target, Vector
        yield from self.selector_range_ref.getInputSchema()
        yield from self.selector_range_target.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return (
            (self.name_mask_ref, Vector),
            (self.name_mask_target, Vector),
            (self.name_count, Scalar),
            (self.name_count_ref, Scalar),
            (self.name_count_target, Scalar),
            (self.name_completeness, Scalar),
            (self.name_completeness_bad_match, Scalar),
            (self.name_completeness_good_match, Scalar),
            (self.name_purity, Scalar),
            (self.name_purity_bad_match, Scalar),
            (self.name_purity_good_match, Scalar),
            (self.name_range_maximum, Scalar),
            (self.name_range_minimum, Scalar),
        )

    def getFormattedOutputKeys(self, **kwargs) -> dict[str, str]:
        """Return the mapping from unformatted output schema keys to formatted.

        Returns
        -------
        result : dict[`str`, `str`]
            A dict with formatted key values for unformatted keys.
        """
        has_band = "band" in kwargs
        kwargs_format = {}
        if has_band:
            kwargs_format["band"] = kwargs["band"]

        return {
            self.name_mask_ref: self.name_mask_ref.format(**kwargs_format),
            self.name_mask_target: self.name_mask_target.format(**kwargs_format),
            self.name_count: self.name_count.format(**kwargs_format),
            self.name_count_ref: self.name_count_ref.format(**kwargs_format),
            self.name_count_target: self.name_count_target.format(**kwargs_format),
            self.name_completeness: self.name_completeness.format(**kwargs_format),
            self.name_completeness_bad_match: self.name_completeness_bad_match.format(**kwargs_format),
            self.name_completeness_good_match: self.name_completeness_good_match.format(**kwargs_format),
            self.name_purity: self.name_purity.format(**kwargs_format),
            self.name_purity_bad_match: self.name_purity_bad_match.format(**kwargs_format),
            self.name_purity_good_match: self.name_purity_good_match.format(**kwargs_format),
            self.name_range_maximum: self.name_range_maximum.format(**kwargs_format),
            self.name_range_minimum: self.name_range_minimum.format(**kwargs_format),
        }

    @property
    def name_count(self):
        return f"{self.name_prefix}count{self.name_suffix}"

    @property
    def name_count_ref(self):
        return f"{self.name_prefix}count_ref{self.name_suffix}"

    @property
    def name_count_target(self):
        return f"{self.name_prefix}count_target{self.name_suffix}"

    @property
    def name_mask_ref(self):
        return f"{self.name_prefix}mask_ref{self.name_suffix}"

    @property
    def name_mask_target(self):
        return f"{self.name_prefix}mask_ref{self.name_suffix}"

    @property
    def name_completeness(self):
        return f"{self.name_prefix}completeness{self.name_suffix}"

    @property
    def name_completeness_bad_match(self):
        return f"{self.name_prefix}completeness_bad_match{self.name_suffix}"

    @property
    def name_completeness_good_match(self):
        return f"{self.name_prefix}completeness_good_match{self.name_suffix}"

    @property
    def name_purity(self):
        return f"{self.name_prefix}purity{self.name_suffix}"

    @property
    def name_purity_bad_match(self):
        return f"{self.name_prefix}purity_bad_match{self.name_suffix}"

    @property
    def name_purity_good_match(self):
        return f"{self.name_prefix}purity_good_match{self.name_suffix}"

    @property
    def name_range_maximum(self):
        return f"{self.name_prefix}range_maximum{self.name_suffix}"

    @property
    def name_range_minimum(self):
        return f"{self.name_prefix}range_minimum{self.name_suffix}"

    def name_mag_completeness(self, name_threshold: str):
        name_threshold = f"_{name_threshold}" if name_threshold else name_threshold
        return f"{self.name_prefix}mag_completeness{name_threshold}{self.name_suffix}"

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        if (self.selector_range_ref.minimum != self.selector_range_target.minimum) or (
            self.selector_range_ref.maximum != self.selector_range_target.maximum
        ):
            raise ValueError(
                f"{self.selector_range_ref.minimum=} != {self.selector_range_target.minimum=} or"
                f" {self.selector_range_ref.maximum=} != {self.selector_range_target.maximum=}"
            )
        results = {}
        mask_kw = kwargs.get("mask")
        mask_ref = self.selector_range_ref(data, **kwargs)
        mask_target = self.selector_range_target(data, **kwargs)
        for mask, key_new in ((mask_ref, self.key_mask_ref), (mask_target, self.key_mask_target)):
            if key_new:
                mask_new = data[key_new]
                if mask_kw:
                    mask_new = mask[mask_kw]
                mask &= mask_new

        results[self.name_mask_ref] = mask_ref
        results[self.name_mask_target] = mask_target

        n_ref = np.sum(mask_ref)
        n_target = np.sum(mask_target)
        mask_any = mask_ref | mask_target
        matched = data[self.key_match_distance] >= 0
        if mask_kw:
            matched = matched[mask_kw]

        matched_ref = matched & mask_ref
        matched_target = matched & mask_target
        n_matched_ref = np.sum(matched_ref)
        n_matched_target = np.sum(matched & mask_target)

        matched_good = data[self.key_matched_class]
        if mask_kw:
            matched_good = matched_good[mask_kw]

        n_matched_same_ref = np.sum(matched_good & matched_ref)
        n_matched_same_target = np.sum(matched_good & matched_target)

        results[self.name_count] = np.sum(mask_any)
        results[self.name_count_ref] = n_ref
        results[self.name_count_target] = n_target
        results[self.name_completeness] = divide(n_matched_ref, n_ref)
        results[self.name_completeness_bad_match] = divide(n_matched_ref - n_matched_same_ref, n_ref)
        results[self.name_completeness_good_match] = divide(n_matched_same_ref, n_ref)
        results[self.name_purity] = divide(n_matched_target, n_target)
        results[self.name_purity_bad_match] = divide(n_matched_target - n_matched_same_target, n_target)
        results[self.name_purity_good_match] = divide(n_matched_same_target, n_target)

        results[self.name_range_maximum] = self.selector_range_ref.maximum
        results[self.name_range_minimum] = self.selector_range_ref.minimum

        keys_formatted = self.getFormattedOutputKeys(**kwargs)
        results = {key_new: results[key_old] for key_old, key_new in keys_formatted.items()}

        return results
