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
    "Completeness",
    "TargetInjectedCatDelta",
    "TargetInjectedCatDeltaScatterAstrom",
    "TargetInjectedCatDeltaRAScatterPlot",
    "TargetInjectedCatDeltaDecScatterPlot",
    "TargetInjectedCatDeltaScatterPhotom",
    "TargetInjectedCatDeltaPsfScatterPlot",
    "TargetInjectedCatDeltaCModelScatterPlot",
    "TargetInjectedCatDeltaSkyPlot",
    "TargetInjectedCatDeltaSkyPlotAstrom",
    "TargetInjectedCatDeltaRASkyPlot",
    "TargetInjectedCatDeltaDecSkyPlot",
    "TargetInjectedCatDeltaSkyPlotPhotom",
    "TargetInjectedCatDeltaPsfSkyPlot",
    "TargetInjectedCatDeltaCModelSkyPlot",
    "TargetInjectedCatDeltaMetrics",
)

from ..actions.plot.completenessPlot import CompletenessHist
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.plot.skyPlot import SkyPlot

from ..actions.scalar.scalarActions import Mag50Action, MedianAction, SigmaMadAction
from ..actions.vector import (
    AngularSeparation,
    ConvertFluxToMag,
    ConvertUnits,
    DownselectVector,
    LoadVector,
    MagDiff,
    RAcosDec,
    RangeSelector,
    SnSelector,
    StarSelector,
    SubtractVector,
)
from ..contexts import CoaddContext, RefMatchContext
from ..interfaces import AnalysisTool
from lsst.pex.config import Field


class Completeness(AnalysisTool):
    """Plot the fraction of injected sources recovered by input magnitude."""

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness_target"

        self.process.buildActions.mag = LoadVector()
        self.process.buildActions.mag.vectorKey = "{band}_mag_ref"
        # self.process.buildActions.mag = ConvertFluxToMag()
        # self.process.buildActions.mag.vectorKey = "{band}_psfFlux_target"
        self.process.buildActions.matchDistance = LoadVector()
        self.process.buildActions.matchDistance.vectorKey = "matchDistance"

        self.process.calculateActions.mag50 = Mag50Action()
        self.process.calculateActions.mag50.vectorKey = "{band}_mag_ref"
        self.process.calculateActions.mag50.matchDistanceKey = "matchDistance"

        self.produce.plot = CompletenessHist()
        self.produce.metric.units = {"mag50": "mag"}
        self.produce.metric.newNames = {"mag50": "{band}_mag50"}

        self.applyContext(RefMatchContext)


class TargetInjectedCatDelta(AnalysisTool):
    """Plot the difference between a target catalog and an
    injected catalog for the quantity set in `setDefaults`.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        """Apply coadd options for the ref cat plots.
        Applies the coadd plot flag selector and sets
        flux types.
        """
        # self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness_target"

        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.lowSNSelector.threshold = 1
        self.process.calculateActions.stars.highSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.highSNSelector.threshold = 100
        self.process.calculateActions.stars.fluxType = "{band}_psfFlux_target"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux_target"
        self.process.buildActions.patch = LoadVector()
        self.process.buildActions.patch.vectorKey = "patch_target"


class TargetInjectedCatDeltaScatterAstrom(TargetInjectedCatDelta):
    """Plot the difference in milliseconds between a target catalog and an
    injected catalog for the coordinate set in `setDefaults`. Plot it on
    a scatter plot.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )
        self.process.buildActions.xStars = ConvertFluxToMag()
        self.process.buildActions.xStars.vectorKey = "{band}_psfFlux_target"
        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.fluxType = "{band}_psfFlux_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetInjectedCatDeltaRAScatterPlot(TargetInjectedCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and an injected catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_target"
        )
        self.process.buildActions.yStars.buildAction.actionB = RAcosDec(raKey="ra_ref", decKey="dec_ref")

        self.produce.yAxisLabel = "RA$_{{output}}$ - RA$_{{input}}$ (marcsec)"


class TargetInjectedCatDeltaDecScatterPlot(TargetInjectedCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and an injected catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.yStars.buildAction.actionB = LoadVector(vectorKey="dec_ref")

        self.produce.yAxisLabel = "Dec$_{{output}}$ - Dec$_{{input}}$ (marcsec)"


class TargetInjectedCatDeltaScatterPhotom(TargetInjectedCatDelta):
    """Plot the difference in millimags between a target catalog and an
    injected catalog for the flux type set in `setDefaults`.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars = MagDiff()
        self.process.buildActions.yStars.col2 = "{band}_mag_ref"
        self.process.buildActions.yStars.fluxUnits2 = "mag(AB)"

        self.process.buildActions.xStars = ConvertFluxToMag()
        self.process.buildActions.xStars.vectorKey = "{band}_psfFlux_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.yAxisLabel = "Output Mag - Input Mag (mmag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetInjectedCatDeltaPsfScatterPlot(TargetInjectedCatDeltaScatterPhotom):
    """Plot the difference in millimags between the PSF flux
    of a target catalog and an injected catalog
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "{band}_psfFlux_target"


class TargetInjectedCatDeltaCModelScatterPlot(TargetInjectedCatDeltaScatterPhotom):
    """Plot the difference in millimags between the CModel flux
    of a target catalog and an injected catalog.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "{band}_cModelFlux_target"


class TargetInjectedCatDeltaSkyPlot(TargetInjectedCatDelta):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference of the quantity defined in the vector key parameter between
    the target and injected catalog as the color.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra_target"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec_target"

        self.produce = SkyPlot()
        self.produce.plotTypes = ["stars"]
        self.produce.xAxisLabel = "R.A. (degrees)"
        self.produce.yAxisLabel = "Dec. (degrees)"
        self.produce.plotOutlines = False


class TargetInjectedCatDeltaSkyPlotAstrom(TargetInjectedCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and injected catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )

        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)

        # self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux_target"


class TargetInjectedCatDeltaRASkyPlot(TargetInjectedCatDeltaSkyPlotAstrom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the RA of the target and injected catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.zStars.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_target"
        )
        self.process.buildActions.zStars.buildAction.actionB = RAcosDec(raKey="ra_ref", decKey="dec_ref")

        self.produce.plotName = "injected_astromDiffSky_RA"
        self.produce.zAxisLabel = "RA$_{{output}}$ - RA$_{{input}}$ (marcsec)"


class TargetInjectedCatDeltaDecSkyPlot(TargetInjectedCatDeltaSkyPlotAstrom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the Dec of the target and injected catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.zStars.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.zStars.buildAction.actionB = LoadVector(vectorKey="dec_ref")

        self.produce.plotName = "injected_astromDiffSky_Dec"
        self.produce.zAxisLabel = "Dec$_{{output}}$ - Dec$_{{input}}$ (marcsec)"


class TargetInjectedCatDeltaSkyPlotPhotom(TargetInjectedCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the photometry of the target and injected catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = MagDiff()
        self.process.buildActions.zStars.col2 = "{band}_mag_ref"
        self.process.buildActions.zStars.fluxUnits2 = "mag(AB)"

        self.produce.plotName = "injected_photomDiffSky_{band}"
        self.produce.zAxisLabel = "Output Mag - Input Mag (mmag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetInjectedCatDeltaPsfSkyPlot(TargetInjectedCatDeltaSkyPlotPhotom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the PSF photometry of the target and injected
    catalog as the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars.col1 = "{band}_psfFlux_target"


class TargetInjectedCatDeltaCModelSkyPlot(TargetInjectedCatDeltaSkyPlotPhotom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the CModel photometry of the target and injected
    catalog as the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars.col1 = "{band}_cModelFlux_target"


class TargetInjectedCatDeltaMetrics(AnalysisTool):
    """Calculate the AA1 metric and the sigma MAD from the difference between
    the target and injected catalog coordinates.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        """Apply coadd options for the metrics. Applies the coadd plot flag
        selector and sets flux types.
        """
        # self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness_target"

        self.applyContext(RefMatchContext)

        self.process.buildActions.mags.vectorKey = "{band}_psfFlux_target"

        self.produce.metric.newNames = {
            "injected_AA1_RA": "injected_{band}_AA1_RA_coadd",
            "injected_AA1_sigmaMad_RA": "injected_{band}_AA1_sigmaMad_RA_coadd",
            "injected_AA1_Dec": "injected_{band}_AA1_Dec_coadd",
            "injected_AA1_sigmaMad_Dec": "injected_{band}_AA1_sigmaMad_Dec_coadd",
            "injected_AA1_tot": "injected_{band}_AA1_tot_coadd",
            "injected_AA1_sigmaMad_tot": "injected_{band}_AA1_sigmaMad_tot_coadd",
        }

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.starSelector = StarSelector()

        # Calculate difference in RA
        self.process.buildActions.astromDiffRA = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )
        self.process.buildActions.astromDiffRA.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="dec_ref"
        )
        self.process.buildActions.astromDiffRA.buildAction.actionB = RAcosDec(
            raKey="ra_ref", decKey="dec_ref"
        )
        # Calculate difference in Dec
        self.process.buildActions.astromDiffDec = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )
        self.process.buildActions.astromDiffDec.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.astromDiffDec.buildAction.actionB = LoadVector(vectorKey="dec_ref")
        # Calculate total difference (using astropy)
        self.process.buildActions.astromDiffTot = AngularSeparation(
            raKey_A="coord_ra_target",
            decKey_A="coord_dec_target",
            raKey_B="ra_ref",
            decKey_B="dec_ref",
        )

        self.process.buildActions.mags = ConvertFluxToMag()

        # Filter down to only objects with mag 17-21.5
        self.process.filterActions.brightStarsRA = DownselectVector(vectorKey="astromDiffRA")
        self.process.filterActions.brightStarsRA.selector = RangeSelector(
            vectorKey="mags", minimum=17, maximum=21.5
        )

        self.process.filterActions.brightStarsDec = DownselectVector(vectorKey="astromDiffDec")
        self.process.filterActions.brightStarsDec.selector = RangeSelector(
            vectorKey="mags", minimum=17, maximum=21.5
        )

        self.process.filterActions.brightStarsTot = DownselectVector(vectorKey="astromDiffTot")
        self.process.filterActions.brightStarsTot.selector = RangeSelector(
            vectorKey="mags", minimum=17, maximum=21.5
        )

        # Calculate median and sigmaMad
        self.process.calculateActions.injected_AA1_RA = MedianAction(vectorKey="brightStarsRA")
        self.process.calculateActions.injected_AA1_sigmaMad_RA = SigmaMadAction(vectorKey="brightStarsRA")

        self.process.calculateActions.injected_AA1_Dec = MedianAction(vectorKey="brightStarsDec")
        self.process.calculateActions.injected_AA1_sigmaMad_Dec = SigmaMadAction(vectorKey="brightStarsDec")

        self.process.calculateActions.injected_AA1_tot = MedianAction(vectorKey="brightStarsTot")
        self.process.calculateActions.injected_AA1_sigmaMad_tot = SigmaMadAction(vectorKey="brightStarsTot")

        self.produce.metric.units = {
            "injected_AA1_RA": "mas",
            "injected_AA1_sigmaMad_RA": "mas",
            "injected_AA1_Dec": "mas",
            "injected_AA1_sigmaMad_Dec": "mas",
            "injected_AA1_tot": "mas",
            "injected_AA1_sigmaMad_tot": "mas",
        }
