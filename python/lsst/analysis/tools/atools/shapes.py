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
    "BasePsfResidualCoadd",
    "BasePsfResidualVisit",
    "E1DiffScatter",
    "E1DiffScatterVisit",
    "E2DiffScatter",
    "E2DiffScatterVisit",
    "ShapeSizeFractionalDiffScatter",
    "ShapeSizeFractionalDiffScatterVisit",
    "E1DiffSky",
    "E1DiffSkyVisit",
    "E2DiffSky",
    "E2DiffSkyVisit",
    "ShapeSizeFractionalDiffSky",
    "ShapeSizeFractionalDiffSkyVisit",
    "RhoStatistics",
    "RelativeSizeResidualPlot",
    "EBase",
    "EScatter",
    "E1ScatterVisit",
    "E2ScatterVisit",
    "ESky",
    "E1SkyVisit",
    "E2SkyVisit",
    "EFocalPlane",
    "E1FocalPlane",
    "E2FocalPlane",
    "ShapeSizeFractionalDiffFocalPlane",
)

from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ..actions.keyedData import KeyedScalars
from ..actions.plot import FocalPlanePlot
from ..actions.plot.rhoStatisticsPlot import RhoStatisticsPlot
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.plot.skyPlot import SkyPlot
from ..actions.scalar import CountAction, MedianAction, SigmaMadAction
from ..actions.vector import (
    CalcE,
    CalcE1,
    CalcE2,
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
    VisitPlotFlagSelector,
)
from ..contexts import CoaddContext, VisitContext
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


class BasePsfResidualCoadd(AnalysisTool):
    """Shared configuration for `prep` and `process` stages of PSF residuals.

    This is a mixin class used by `BasePsfResidualScatterPlot` and
    `BasePsfResidualMetric` to share common default configuration.
    """

    parameterizedBand: bool = True  # False

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)

        self.process.calculateActions.stars = ScatterPlotStatsAction(
            vectorKey="yStars",
        )
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

        self.process.buildActions.patch = LoadVector()
        self.process.buildActions.patch.vectorKey = "patch"

        self.process.buildActions.xStars = ConvertFluxToMag(vectorKey="{band}_psfFlux")
        self.process.buildActions.psfFlux = LoadVector(vectorKey="{band}_psfFlux")
        self.process.buildActions.psfFluxErr = LoadVector(vectorKey="{band}_psfFluxErr")


class BasePsfResidualVisit(AnalysisTool):
    """Shared configuration for `prep` and `process` stages of PSF residuals.

    This is a mixin class used by `BasePsfResidualScatterPlot` and
    `BasePsfResidualMetric` to share common default configuration.
    """

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.flagSelector = VisitPlotFlagSelector()
        self.prep.selectors.snSelector = SnSelector(fluxType="psfFlux", threshold=50)
        self.prep.selectors.starSelector = StarSelector(vectorKey="extendedness")

        self.process.buildActions.xStars = ConvertFluxToMag(vectorKey="psfFlux")

        self.process.buildActions.mags = ConvertFluxToMag(vectorKey="psfFlux")

        self.process.buildActions.fracDiff = FractionalDifference(
            actionA=CalcMomentSize(colXx="ixx", colYy="iyy", colXy="ixy"),
            actionB=CalcMomentSize(colXx="ixxPSF", colYy="iyyPSF", colXy="ixyPSF"),
        )
        # Define an eDiff action and let e1Diff and e2Diff differ only in
        # component.
        self.process.buildActions.eDiff = CalcEDiff(
            colA=CalcE(colXx="ixx", colYy="iyy", colXy="ixy"),
            colB=CalcE(colXx="ixxPSF", colYy="iyyPSF", colXy="ixyPSF"),
        )
        self.process.buildActions.e1Diff = self.process.buildActions.eDiff
        self.process.buildActions.e1Diff.component = "1"
        self.process.buildActions.e2Diff = self.process.buildActions.eDiff
        self.process.buildActions.e2Diff.component = "2"


class SkyCoadd(BasePsfResidualCoadd):
    """A base class for coadd level sky plots."""

    parameterizedBand: bool = True  # False

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec"

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"

        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["stars"]
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"


class SkyVisit(BasePsfResidualVisit):
    """A base class for visit level sky plots."""

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec"

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "psfFlux"

        self.produce.plot = SkyPlot()

        self.produce.plot.plotTypes = ["stars"]
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"

        self.produce.plot.plotOutlines = False


class ScatterCoadd(BasePsfResidualCoadd):
    """A base class for coadd level scatter plots."""

    parameterizedBand: bool = True  # False

    def setDefaults(self):
        super().setDefaults()
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


class ScatterVisit(BasePsfResidualVisit):
    """A base class for visit level sky plots."""

    def setDefaults(self):
        super().setDefaults()
        self.process.calculateActions.stars = ScatterPlotStatsAction(
            vectorKey="yStars",
        )

        # use the downselected psfFlux
        self.process.calculateActions.stars.highSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.highSNSelector.threshold = 500
        self.process.calculateActions.stars.lowSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.lowSNSelector.threshold = 100
        self.process.calculateActions.stars.fluxType = "psfFlux"

        self.produce.metric.units = {
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        }

        self.produce.metric.newNames = {
            "{band}_highSNStars_median": "highSNStars_median",
            "{band}_highSNStars_sigmaMad": "highSNStars_sigmaMad",
            "{band}_highSNStars_count": "highSNStars_count",
            "{band}_lowSNStars_median": "lowSNStars_median",
            "{band}_lowSNStars_sigmaMad": "lowSNStars_sigmaMad",
            "{band}_lowSNStars_count": "lowSNStars_count",
        }

        self.produce.plot = ScatterPlotWithTwoHists()

        self.produce.plot.plotTypes = ["stars"]
        self.produce.plot.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.plot.magLabel = "PSF Magnitude (mag)"

        self.produce.plot.addSummaryPlot = False


class E1DiffScatter(ScatterCoadd):
    """The difference between E1 and E1 PSF plotted on
    a scatter plot. Used in coaddQualityCore.
    """

    def setDefaults(self):
        super().setDefaults()
        self.applyContext(CoaddContext)
        self.process.filterActions.yStars = LoadVector(vectorKey="e1Diff")
        self.produce.plot.yAxisLabel = "Ellipticity residuals (e1 - e1_PSF)"


class E1DiffScatterVisit(ScatterVisit):
    """The difference between E1 and E1 PSF plotted on
    a scatter plot for visit level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = LoadVector(vectorKey="e1Diff")
        self.applyContext(VisitContext)
        self.produce.plot.yAxisLabel = "Ellipticity residuals (e1 - e1_PSF)"


class E2DiffScatter(ScatterCoadd):
    """The difference between E2 and E2 PSF plotted on
    a scatter plot for coadd level data.
    Used in coaddQualityCore.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = LoadVector(vectorKey="e2Diff")
        self.produce.plot.yAxisLabel = "Ellipticity residuals (e2 - e2_PSF)"
        self.applyContext(CoaddContext)


class E2DiffScatterVisit(ScatterVisit):
    """The difference between E2 and E2 PSF plotted on
    a scatter plot for visit level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = LoadVector(vectorKey="e2Diff")
        self.applyContext(VisitContext)
        self.produce.plot.yAxisLabel = "Ellipticity residuals (e2 - e2_PSF)"


class ShapeSizeFractionalDiffScatter(ScatterCoadd):
    """The difference between shape and shape PSF plotted on
    a scatter plot for coadd level data.
    Used in coaddQualityCore.
    """

    def setDefaults(self):
        super().setDefaults()
        self.applyContext(CoaddContext)
        self.process.filterActions.yStars = LoadVector(vectorKey="fracDiff")
        self.produce.plot.yAxisLabel = "Fractional size residuals (S/S_PSF - 1)"


class ShapeSizeFractionalDiffScatterVisit(ScatterVisit):
    """The difference between shape and shape PSF plotted on
    a scatter plot for visit level data.
    Used in debugPsf.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = LoadVector(vectorKey="fracDiff")
        self.applyContext(VisitContext)
        self.produce.plot.yAxisLabel = "Fractional size residuals (S/S_PSF - 1)"


class E1DiffSky(SkyCoadd):
    """The difference between E1 and E1 PSF plotted on
    a sky plot for coadd level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.zStars = LoadVector(vectorKey="e1Diff")
        self.produce.plot.zAxisLabel = "E1 Diff"
        self.produce.plot.plotName = "E1 Diff Sky"
        self.applyContext(CoaddContext)


class E1DiffSkyVisit(SkyVisit):
    """The difference between E1 and E1 PSF plotted on
    a sky plot for visit level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.zStars = LoadVector(vectorKey="e1Diff")
        self.produce.plot.zAxisLabel = "E1 Diff"
        self.produce.plot.plotName = "E1 Diff Sky"
        self.applyContext(VisitContext)


class E2DiffSky(SkyCoadd):
    """The difference between E2 and E2 PSF plotted on
    a sky plot for coadd level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.zStars = LoadVector(vectorKey="e2Diff")
        self.produce.plot.zAxisLabel = "E2 Diff"
        self.produce.plot.plotName = "E2 Diff Sky"
        self.applyContext(CoaddContext)


class E2DiffSkyVisit(SkyVisit):
    """The difference between E2 and E2 PSF plotted on
    a sky plot for visit level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.zStars = LoadVector(vectorKey="e2Diff")
        self.produce.plot.zAxisLabel = "E2 Diff"
        self.produce.plot.plotName = "E2 Diff Sky"
        self.applyContext(VisitContext)


class ShapeSizeFractionalDiffSky(SkyCoadd):
    """The difference between shape and shape PSF plotted on
    a sky plot for coadd level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.applyContext(CoaddContext)
        self.process.filterActions.zStars = LoadVector(vectorKey="fracDiff")
        self.produce.plot.zAxisLabel = "E1 Diff"
        self.produce.plot.plotName = "Fractional Diff Sky"


class ShapeSizeFractionalDiffSkyVisit(SkyVisit):
    """The difference between shape and shape PSF plotted on
    a sky plot for visit level data.
    Used in debugPSF.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.zStars = LoadVector(vectorKey="fracDiff")
        self.produce.plot.zAxisLabel = "E1 Diff"
        self.produce.plot.plotName = "Fractional Diff Sky"

        self.applyContext(VisitContext)


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


class RelativeSizeResidualPlot(AnalysisTool):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()
        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "extendedness"
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "x"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "y"
        self.process.buildActions.z = FractionalDifference(
            actionA=CalcMomentSize(colXx="ixx", colYy="iyy", colXy="ixy", sizeType="trace"),
            actionB=CalcMomentSize(colXx="ixxPSF", colYy="iyyPSF", colXy="ixyPSF", sizeType="trace"),
        )
        self.process.buildActions.detector = LoadVector(vectorKey="detector")
        self.produce.plot = FocalPlanePlot()
        self.produce.plot.zAxisLabel = "Residuals"
        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.threshold = 20
        self.process.buildActions.statMask.fluxType = "psfFlux"


class EBase(AnalysisTool):
    """A base class for plotting ellipticity plots."""

    def setDafaults(self):
        super().setDefaults()

        self.prep.selectors.snSelector = SnSelector()

        self.prep.selectors.starSelector = StarSelector()

        self.produce.plot.plotTypes = ["stars"]

        def visitContext(self):
            self.parameterizedBand = False
            self.prep.selectors.flagSelector = VisitPlotFlagSelector()

            self.prep.selectors.snSelector.fluxType = "psfFlux"
            self.prep.selectors.snSelector.threshold = 100

            self.prep.selectors.starSelector.vectorKey = "extendedness"

            self.process.buildActions.xStars.vectorKey = "psfFlux"


class EScatter(EBase):
    """Base class for ellipticity scatter plots."""

    def setDafaults(self):
        super().setDefaults()

        self.process.buildActions.xStars = ConvertFluxToMag()

        self.process.calculateActions.median = MedianAction()
        self.process.calculateActions.median.vectorKey = "yStars"
        self.process.calculateActions.sigmaMad = SigmaMadAction()
        self.process.calculateActions.sigmaMad.vectorKey = "yStars"

        self.produce.plot = ScatterPlotWithTwoHists()
        self.produce.plot.xAxisLabel = "PSF Magnitude (Mag)"
        self.produce.plot.magLabel = "PSF Magnitude (Mag)"

        self.produce.metric.units = {"median": "pix", "sigmaMad": "pix"}


class E1ScatterVisit(EScatter):
    """Visit level scatter plot for E1."""

    parameterizedBand: bool = False

    def setDafaults(self):
        super().setDefaults()

        self.process.buildActions.yStarsAll = CalcE1(colXx="ixx", colYy="iyy", colXy="ixy")
        self.process.buildActions.yStars = DownselectVector(
            vectorKey="yStarsAll", selector=VectorSelector(vectorKey="starSelector")
        )

        self.produce.plot.yAxisLabel = "e1"

        self.produce.metric.newNames = {
            "median": "e1_median",
            "sigmaMad": "e1_sigmaMad",
        }
        self.applyContext(VisitContext)


class E2ScatterVisit(EScatter):
    """Visit level scatter plot for E2."""

    parameterizedBand: bool = False

    def setDafaults(self):
        super().setDefaults()

        self.process.buildActions.yStarsAll = CalcE2(colXx="ixx", colYy="iyy", colXy="ixy")
        self.process.buildActions.yStars = DownselectVector(
            vectorKey="yStarsAll", selector=VectorSelector(vectorKey="starSelector")
        )

        self.produce.plot.yAxisLabel = "e2"

        self.produce.metric.newNames = {
            "median": "e2_median",
            "sigmaMad": "e2_sigmaMad",
        }
        self.applyContext(VisitContext)


class ESky(EBase):
    """Base class for ellipticity sky plots."""

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec"

        self.produce.plot = SkyPlot()
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"


class E1SkyVisit(ESky):
    """A sky plot of E1 for visit level data."""

    parameterizedBand: bool = False

    def setDafaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = CalcE1(colXx="ixx", colYy="iyy", colXy="ixy")

        self.produce.plot.zAxisLabel = "e1"

        self.applyContext(VisitContext)


class E2SkyVisit(ESky):
    """A visit level sky plot for E2."""

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = CalcE2(colXx="ixx", colYy="iyy", colXy="ixy")

        self.produce.plot.zAxisLabel = "e2"

        self.applyContext(VisitContext)


class EFocalPlane(EBase):
    """Base class for ellipticity focal plane plots."""

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.x = LoadVector(vectorKey="x")
        self.process.buildActions.y = LoadVector(vectorKey="y")
        self.process.buildActions.detector = LoadVector(vectorKey="detector")

        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.threshold = 20
        self.process.buildActions.statMask.fluxType = "psfFlux"

        self.produce.plot = FocalPlanePlot()


class E1FocalPlane(EFocalPlane):
    """A focal plane plot of E1."""

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.z = CalcE1(colXx="ixx", colYy="iyy", colXy="ixy")
        self.produce.plot.zAxisLabel = "e1"


class E2FocalPlane(EFocalPlane):
    """A focal plane plot of E2."""

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.z = CalcE2(colXx="ixx", colYy="iyy", colXy="ixy")
        self.produce.plot.zAxisLabel = "e2"


class ShapeSizeFractionalDiffFocalPlane(EFocalPlane):
    """A focal plane plot of shape - psf shape."""

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.z = FractionalDifference(
            actionA=CalcMomentSize(colXx="ixx", colYy="iyy", colXy="ixy"),
            actionB=CalcMomentSize(colXx="ixxPSF", colYy="iyyPSF", colXy="ixyPSF"),
        )
        self.produce.plot.zAxisLabel = "Fractional size residuals (S/S_PSF - 1)"
