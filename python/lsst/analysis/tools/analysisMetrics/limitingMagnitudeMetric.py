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

__all__ = ("FiveSigmaPointSourceDepthMetric",)

from ..actions.scalar.scalarActions import MeanAction, MedianAction
from ..actions.vector.selectors import SnSelector, StarSelector
from ..actions.vector.vectorActions import MagColumnNanoJansky
from ..interfaces import AnalysisMetric


class FiveSigmaPointSourceDepthMetric(AnalysisMetric):
    """Calculate the five-sigma point source depth of a visit, based on the
    PSF flux and its reported error. By default the calculation selects
    objects between 4.75 < S/N < 5.25, but these limits are configurable.
    The flux type to use for selection is also configurable. Both the median
    and mean 5-sigma depths are returned.
    """

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "extendedness"

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        # Select between 4.75 < SNR < 5.25 to get a decent sample:
        self.prep.selectors.snSelector.threshold = 4.75
        self.prep.selectors.snSelector.maxSN = 5.25

        self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="psfFlux")
        self.process.calculateActions.median5sigmaDepth = MedianAction(vectorKey="mags")
        self.process.calculateActions.mean5sigmaDepth = MeanAction(vectorKey="mags")

        self.produce.units = {  # type: ignore
            "median5sigmaDepth": "mag",
            "mean5sigmaDepth": "mag",
        }
