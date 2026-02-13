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

__all__ = ("CalcCompletenessHistogramAction", "MagnitudeCompletenessConfig")

import copy

import numpy as np

import lsst.pex.config as pexConfig
from lsst.pex.config.configurableActions import ConfigurableActionField

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema
from ...math import isPercent
from ..config import MagnitudeBinConfig
from .calcBinnedCompleteness import CalcBinnedCompletenessAction


class MagnitudeCompletenessConfig(pexConfig.Config):
    """Configuration for measuring magnitudes at given completeness
    thresholds."""

    completeness_mag_max = pexConfig.Field[float](
        doc="Brightest magnitude to consider checking if completeness is below a percentile threshold for",
        default=18,
    )
    completeness_percentiles = pexConfig.ListField[float](
        doc="The percentiles to find the magnitude at.",
        default=[90.0, 80.0, 50.0],
        itemCheck=isPercent,
    )


class CalcCompletenessHistogramAction(KeyedDataAction):
    """Action to calculate a histogram of completeness vs magnitude."""

    action = ConfigurableActionField[CalcBinnedCompletenessAction](
        doc="The action to compute completeness/purity",
    )
    bins = pexConfig.ConfigField[MagnitudeBinConfig](
        doc="The magnitude bin configuration",
    )
    config_metrics = pexConfig.ConfigField[MagnitudeCompletenessConfig](doc="Metric definition configuration")

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        band = kwargs.get("band")
        bins = tuple(x / 1000.0 for x in reversed(self.bins.get_bins()))
        bin_width = self.bins.mag_width / 1000.0
        n_bins = len(bins)
        keys_raw = {
            key_formatted: key
            for key, key_formatted in self.action.getFormattedOutputKeys(**kwargs).items()
            if not self._is_action_key_mask(key)
        }
        results = {key_formatted: np.zeros(n_bins, dtype=float) for key_formatted in keys_raw.keys()}
        action = copy.copy(self.action)
        name_completeness = action.name_completeness
        if band is not None:
            name_completeness = name_completeness.format(band=band)

        percentile_completeness = {pc: np.nan for pc in sorted(self.config_metrics.completeness_percentiles)}
        percentiles = list(percentile_completeness.keys())
        n_percentiles = len(percentiles)
        idx_percentile = 0
        # isfinite won't work on None
        completeness_last, median_last = np.nan, np.nan

        for idx_rev, minimum in enumerate(bins):
            maximum = minimum + bin_width
            median = (maximum + minimum) / 2.0
            action.selector_range_ref.minimum = minimum
            action.selector_range_ref.maximum = maximum
            action.selector_range_target.minimum = minimum
            action.selector_range_target.maximum = maximum
            result = action(data, **kwargs)
            for key_formatted, array in results.items():
                value = result[key_formatted]
                # The implicit float conversion will generate a warning if the
                # value is masked, so check for that first
                array[n_bins - idx_rev - 1] = np.nan if np.ma.is_masked(value) else value
            if median >= self.config_metrics.completeness_mag_max:
                completeness = result[name_completeness] * 100
                if completeness:
                    if idx_percentile > 0:
                        if completeness < percentiles[idx_percentile - 1]:
                            idx_percentile -= 1
                    while (idx_percentile < n_percentiles) and (completeness > percentiles[idx_percentile]):
                        # Crude linear interpolation
                        # TODO: Replace with a spline or anything better
                        if np.isfinite(completeness_last) and np.isfinite(median_last):
                            percentile = percentiles[idx_percentile]
                            width = completeness - completeness_last
                            # The abs should be unnecessary but just in case
                            magnitude = (
                                abs(completeness - percentile) / width * median_last
                                + abs(percentile - completeness_last) / width * median
                            )
                        else:
                            magnitude = median
                        percentile_completeness[percentiles[idx_percentile]] = magnitude
                        idx_percentile += 1
                    completeness_last, median_last = completeness, median

        for percentile, magnitude in percentile_completeness.items():
            name_percentile = self.getPercentileName(percentile)
            key = action.name_mag_completeness(name_percentile)
            if band is not None:
                key = key.format(band=band)
            results[key] = magnitude

        return results

    def _is_action_key_mask(self, key: str):
        is_mask = key.startswith(f"{self.action.name_prefix}mask")
        return is_mask

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        result = {
            (key, typ) for key, typ in self.action.getOutputSchema() if not self._is_action_key_mask(key)
        }
        return result

    def getPercentileName(self, percentile: float) -> str:
        return f"{percentile:.2f}_pct".replace(".", "p")
