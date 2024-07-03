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

__all__ = ("MagPercentileAction",)

import logging

import numpy as np
from astropy import units as u
from lsst.pex.config import Field, ListField

from ...interfaces import KeyedData, KeyedDataSchema, Scalar, Vector, VectorAction
from ...math import fluxToMag, isPercent

_LOG = logging.getLogger(__name__)


class MagPercentileAction(VectorAction):
    """Calculates the magnitude at the given percentile for completeness"""

    matchDistanceKey = Field[str]("Match distance Vector")
    vectorKey = Field[str](doc="Key of vector which should be loaded")
    fluxUnits = Field[str](doc="Units for the column.", default="nanojansky")
    percentiles = ListField[float](
        doc="The percentiles to find the magnitude at.", default=[16.0, 50.0, 84.0], itemCheck=isPercent
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return (
            (self.matchDistanceKey, Vector),
            (self.vectorKey, Vector),
        )

    def getOutputSchema(self) -> KeyedDataSchema:
        result = []
        for pct in self.percentiles:
            name = self.getPercentileName(pct)
            result.append((name, Scalar))
        return result

    def getPercentileName(self, percentile: float) -> str:
        return f"mag_{percentile:.2f}"

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        matched = np.isfinite(data[self.matchDistanceKey])
        fluxValues = data[self.vectorKey.format(**kwargs)]
        values = fluxToMag(fluxValues, flux_unit=u.Unit(self.fluxUnits))
        nInput, bins = np.histogram(
            values,
            range=(np.nanmin(values), np.nanmax(values)),
            bins=100,
        )
        nOutput, _ = np.histogram(
            values[matched],
            range=(np.nanmin(values[matched]), np.nanmax(values[matched])),
            bins=bins,
        )
        # Find bin where the fraction recovered first falls below a percentile.
        mags: KeyedData = {}
        for pct in self.percentiles:
            name = self.getPercentileName(pct)
            belowPercentile = np.where((nOutput / nInput < pct / 100))[0]
            if len(belowPercentile) == 0:
                mags[name] = np.nan
            else:
                mags[name] = np.min(bins[belowPercentile])
        return mags
