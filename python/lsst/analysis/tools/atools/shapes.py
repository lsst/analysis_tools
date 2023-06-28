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

__all__ = (
    "ShapeSizeFractionalScalars",
    "BasePsfResidual",
    "ShapeSizeFractionalDiff",
    "E1Diff",
    "E2Diff",
    "RhoStatistics",
)

from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ..actions.keyedData import KeyedScalars
from ..actions.plot.rhoStatisticsPlot import RhoStatisticsPlot
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.scalar import CountAction, MedianAction, SigmaMadAction
from ..actions.vector import (
    CalcE,
    CalcEDiff,
    CalcMomentSize,
    CalcRhoStatistics,
    CoaddPlotFlagSelector,
    ConvertFluxToMag,
    DownselectVector,
    FractionalDifference,
    LoadVector,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from ..interfaces import AnalysisTool, KeyedData, VectorAction


class ShapeSizeFractionalScalars(KeyedScalars):
    vectorKey = Field[str](doc="Column key to compute scalars")

    snFluxType = Field[str](doc="column key for the flux type used in SN selection")

    selector = ConfigurableActionField[VectorAction](doc="Selector to use before computing Scalars")

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.sigmaMad = SigmaMadAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.count = CountAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.approxMag = MedianAction(vectorKey=self.snFluxType)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        selection = self.selector(data, **kwargs)
        if mask is not None:
            mask &= selection
        else:
            mask = selection
        return super().__call__(data, **kwargs | dict(mask=mask))


class BasePsfResidual(AnalysisTool):
    """Shared configuration for `prep` and `process` stages of PSF residuals.

    This is a mixin class used by `BasePsfResidualScatterPlot` and
    `BasePsfResidualMetric` to share common default configuration.
    """

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)

        self.process.buildActions.patchWhole = LoadVector()
        self.process.buildActions.patchWhole.vectorKey = "patch"

        self.process.buildActions.mags = ConvertFluxToMag(vectorKey="{band}_psfFlux")

        self.process.buildActions.fracDiff = FractionalDifference(
            actionA=CalcMomentSize(colXx="{band}_ixx", colYy="{band}_iyy", colXy="{band}_ixy"),
            actionB=CalcMomentSize(colXx="{band}_ixxPSF", colYy="{band}_iyyPSF", colXy="{band}_ixyPSF"),
        )
        # Define an eDiff action and let e1Diff and e2Diff differ only in
        # component.
        self.process.buildActions.eDiff = CalcEDiff(
            colA=CalcE(colXx="{band}_ixx", colYy="{band}_iyy", colXy="{band}_ixy"),
            colB=CalcE(colXx="{band}_ixxPSF", colYy="{band}_iyyPSF", colXy="{band}_ixyPSF"),
        )
        self.process.buildActions.e1Diff = self.process.buildActions.eDiff
        self.process.buildActions.e1Diff.component = "1"
        self.process.buildActions.e2Diff = self.process.buildActions.eDiff
        self.process.buildActions.e2Diff.component = "2"
        # pre-compute a stellar selector mask so it can be used in the filter
        # actions while only being computed once, alternatively the stellar
        # selector could be calculated and applied twice in the filter stage
        self.process.buildActions.starSelector = StarSelector()

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

        self.process.calculateActions.stars = ScatterPlotStatsAction(
            vectorKey="yStars",
        )
        # use the downselected psfFlux
        self.process.calculateActions.stars.highSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.lowSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.fluxType = "psfFlux"

        self.produce.metric.units = {
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        }

        self.produce.plot = ScatterPlotWithTwoHists()

        self.produce.plot.plotTypes = ["stars"]
        self.produce.plot.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.plot.magLabel = "PSF Magnitude (mag)"
        self.produce.plot.addSummaryPlot = True


class ShapeSizeFractionalDiff(BasePsfResidual):
    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="fracDiff", selector=VectorSelector(vectorKey="starSelector")
        )
        self.produce.plot.yAxisLabel = "Fractional size residuals (S/S_PSF - 1)"


class E1Diff(BasePsfResidual):
    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="e1Diff", selector=VectorSelector(vectorKey="starSelector")
        )
        self.produce.plot.yAxisLabel = "Ellipticty residuals (e1 - e1_PSF)"


class E2Diff(BasePsfResidual):
    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="e2Diff", selector=VectorSelector(vectorKey="starSelector")
        )
        self.produce.plot.yAxisLabel = "Ellipticty residuals (e2 - e2_PSF)"


class RhoStatistics(AnalysisTool):
    parameterizedBand: bool = True

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)
        self.prep.selectors.starSelector = StarSelector()

        self.process.calculateActions.rho = CalcRhoStatistics()
        self.process.calculateActions.rho.treecorr.nbins = 21
        self.process.calculateActions.rho.treecorr.min_sep = 0.01
        self.process.calculateActions.rho.treecorr.max_sep = 100.0
        self.process.calculateActions.rho.treecorr.sep_units = "arcmin"
        self.process.calculateActions.rho.treecorr.metric = "Arc"

        self.produce.plot = RhoStatisticsPlot()
