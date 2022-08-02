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
    "E1DiffScatterPlot",
    "E2DiffScatterPlot",
    "ShapeSizeFractionalDiffScatterPlot",
    "WPerpPSFPlot",
    "Ap12PsfSkyPlot",
)

from ..actions.keyedData.stellarLocusFit import StellarLocusFitAction
from ..actions.plot.colorColorFitPlot import ColorColorFitPlot
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.plot.skyPlot import SkyPlot
from ..actions.scalar import ApproxFloor
from ..actions.vector import (
    CoaddPlotFlagSelector,
    DownselectVector,
    ExtinctionCorrectedMagDiff,
    LoadVector,
    MagColumnNanoJansky,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from ..analysisParts.shapeSizeFractional import BasePsfResidualMixin
from ..interfaces import AnalysisPlot


class BasePsfResidualScatterPlot(AnalysisPlot, BasePsfResidualMixin):
    """Base class for scatter plots of PSF residuals.

    This is shared by size and ellipticity plots.
    """

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

        self.produce = ScatterPlotWithTwoHists()

        self.produce.plotTypes = ["stars"]
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.addSummaryPlot = True


class ShapeSizeFractionalDiffScatterPlot(BasePsfResidualScatterPlot):
    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="fracDiff", selector=VectorSelector(vectorKey="starSelector")
        )
        self.produce.yAxisLabel = "Fractional size residuals (S/S_PSF - 1)"


class E1DiffScatterPlot(BasePsfResidualScatterPlot):
    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="e1Diff", selector=VectorSelector(vectorKey="starSelector")
        )
        self.produce.yAxisLabel = "Ellipticty residuals (e1 - e1_PSF)"


class E2DiffScatterPlot(BasePsfResidualScatterPlot):
    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="e2Diff", selector=VectorSelector(vectorKey="starSelector")
        )
        self.produce.yAxisLabel = "Ellipticty residuals (e2 - e2_PSF)"


class WPerpPSFPlot(AnalysisPlot):
    # WPerp does not support running in multiband mode
    multiband: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.bands = ["r"]
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 300

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "r_extendedness"

        self.process.buildActions.x = ExtinctionCorrectedMagDiff()
        self.process.buildActions.x.magDiff.col1 = "g_psfFlux"
        self.process.buildActions.x.magDiff.col2 = "r_psfFlux"
        self.process.buildActions.x.magDiff.returnMillimags = False
        self.process.buildActions.y = ExtinctionCorrectedMagDiff()
        self.process.buildActions.y.magDiff.col1 = "r_psfFlux"
        self.process.buildActions.y.magDiff.col2 = "i_psfFlux"
        self.process.buildActions.y.magDiff.returnMillimags = False
        self.process.buildActions.mag = MagColumnNanoJansky(vectorKey="r_psfFlux")

        self.process.calculateActions.approxMagDepth = ApproxFloor(vectorKey="mag")
        self.process.calculateActions.wPerp_psfFlux = StellarLocusFitAction()
        self.process.calculateActions.wPerp_psfFlux.stellarLocusFitDict = {
            "xMin": 0.28,
            "xMax": 1.0,
            "yMin": 0.02,
            "yMax": 0.48,
            "mHW": 0.52,
            "bHW": -0.08,
        }

        self.produce = ColorColorFitPlot()
        self.produce.plotName = "wPerp_psfFlux"


class Ap12PsfSkyPlot(AnalysisPlot):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["{band}"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 300

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness"

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"

        self.process.buildActions.zStars = ExtinctionCorrectedMagDiff()
        self.process.buildActions.zStars.magDiff.col1 = "{band}_ap12Flux"
        self.process.buildActions.zStars.magDiff.col2 = "{band}_psfFlux"

        self.produce = SkyPlot()
        self.produce.plotTypes = ["stars"]
        self.produce.plotName = "ap12-psf_{band}"
        self.produce.xAxisLabel = "R.A. (degrees)"
        self.produce.yAxisLabel = "Dec. (degrees)"
        self.produce.zAxisLabel = "Ap 12 - PSF [mag]"
        self.produce.plotOutlines = False
