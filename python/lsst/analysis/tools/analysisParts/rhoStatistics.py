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

__all__ = ("RhoStatisticsMixin",)

from ..actions.vector import (
    CalcRhoStatistics,
    CoaddPlotFlagSelector,
    DownselectVector,
    FlagSelector,
    LoadVector,
    MagColumnNanoJansky,
    SnSelector,
    VectorSelector,
)
from ..interfaces import AnalysisTool


class RhoStatisticsMixin(AnalysisTool):
    parameterizedBand: bool = True

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)

        self.process.buildActions.patchWhole = LoadVector()
        self.process.buildActions.patchWhole.vectorKey = "patch"

        self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="{band}_psfFlux")
        # pre-compute a stellar selector mask so it can be used in the filter
        # actions while only being computed once, alternatively the stellar
        # selector could be calculated and applied twice in the filter stage
        self.process.buildActions.starSelector = FlagSelector(selectWhenTrue=("{band}_calib_psf_used",))

        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags", selector=VectorSelector(vectorKey="starSelector")
        )
        # downselect the psfFlux as well
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="{band}_psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="{band}_psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.filterActions.patch = DownselectVector(
            vectorKey="patchWhole", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.stars = CalcRhoStatistics()

        self.process.calculateActions.stars.treecorr.nbins = 10
        self.process.calculateActions.stars.treecorr.min_sep = 0.1
        self.process.calculateActions.stars.treecorr.max_sep = 100.0
        self.process.calculateActions.stars.treecorr.sep_units = "arcmin"
