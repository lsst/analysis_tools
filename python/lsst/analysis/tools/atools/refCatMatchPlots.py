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
    "TargetRefCatDeltaRAScatterPlot",
    "TargetRefCatDeltaDecScatterPlot",
    "TargetRefCatDeltaRASkyPlot",
    "TargetRefCatDeltaDecSkyPlot",
    "TargetRefCatDeltaScatterPhotom",
    "TargetRefCatDeltaScatterAstrom",
    "TargetRefCatDeltaSkyPlotPhotom",
    "TargetRefCatDeltaPsfScatterPlot",
    "TargetRefCatDeltaCModelScatterPlot",
    "TargetRefCatDeltaCModelSkyPlot",
    "TargetRefCatDeltaPsfSkyPlot",
    "TargetRefCatDeltaRASkyVisitPlot",
    "TargetRefCatDeltaAp09ScatterVisitPlot",
    "TargetRefCatDeltaPsfScatterVisitPlot",
    "TargetRefCatDeltaAp09SkyVisitPlot",
    "TargetRefCatDeltaPsfSkyVisitPlot",
    "TargetRefCatDeltaDecSkyVisitPlot",
    "TargetRefCatDeltaRAScatterVisitPlot",
    "TargetRefCatDeltaDecScatterVisitPlot",
    "TargetRefCatDeltaMetrics",
)

from lsst.pex.config import Field

from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.plot.skyPlot import SkyPlot
from ..actions.scalar.scalarActions import MedianAction, SigmaMadAction
from ..actions.vector import (
    AngularSeparation,
    CoaddPlotFlagSelector,
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
    VisitPlotFlagSelector,
)
from ..contexts import CoaddContext, RefMatchContext, VisitContext
from ..interfaces import AnalysisTool


class TargetRefCatDelta(AnalysisTool):
    """Plot the difference between a target catalog and a
    reference catalog for the quantity set in `setDefaults`.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        """Apply coadd options for the ref cat plots.
        Applies the coadd plot flag selector and sets
        flux types.
        """
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux_target"
        self.prep.selectors.snSelector.threshold = 200
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness_target"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux_target"
        self.process.buildActions.patch = LoadVector()
        self.process.buildActions.patch.vectorKey = "patch_target"

    def visitContext(self) -> None:
        """Apply visit options for the ref cat plots.
        Applies the visit plot flag selector and sets
        the flux types.
        """
        self.parameterizedBand = False
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()
        self.prep.selectors.starSelector.vectorKey = "extendedness_target"
        self.prep.selectors.snSelector.fluxType = "psfFlux_target"
        self.prep.selectors.snSelector.threshold = 50

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "psfFlux_target"
        self.process.buildActions.starStatMask.threshold = 200

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.starSelector = StarSelector()


class TargetRefCatDeltaScatterAstrom(TargetRefCatDelta):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`. Plot it on
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
        self.process.calculateActions.stars.lowSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.highSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.fluxType = "{band}_psfFlux_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaScatterAstromVisit(TargetRefCatDelta):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`. Plot it on
    a scatter plot.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )
        self.process.buildActions.xStars = ConvertFluxToMag()
        self.process.buildActions.xStars.vectorKey = "psfFlux_target"
        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "psfFlux_target"
        self.process.calculateActions.stars.highSNSelector.fluxType = "psfFlux_target"
        self.process.calculateActions.stars.fluxType = "psfFlux_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.addSummaryPlot = False
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.applyContext(VisitContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaScatterPhotom(TargetRefCatDelta):
    """Plot the difference in millimags between a target catalog and a
    reference catalog for the flux type set in `setDefaults`.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars = MagDiff()
        self.process.buildActions.yStars.col2 = "{band}_mag_ref"
        self.process.buildActions.yStars.fluxUnits2 = "mag(AB)"

        self.process.buildActions.xStars = ConvertFluxToMag()
        self.process.buildActions.xStars.vectorKey = "{band}_psfFlux_target"

        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.highSNSelector.fluxType = "{band}_psfFlux_target"
        self.process.calculateActions.stars.fluxType = "{band}_psfFlux_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.yAxisLabel = "Output Mag - Ref Mag (mmag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaScatterPhotomVisit(TargetRefCatDelta):
    """Plot the difference in millimags between a target catalog and a
    reference catalog for the flux type set in `setDefaults`.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars = MagDiff()
        self.process.buildActions.yStars.col2 = "mag_ref"
        self.process.buildActions.yStars.fluxUnits2 = "mag(AB)"

        self.process.buildActions.xStars = ConvertFluxToMag()
        self.process.buildActions.xStars.vectorKey = "psfFlux_target"

        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "psfFlux_target"
        self.process.calculateActions.stars.highSNSelector.fluxType = "psfFlux_target"
        self.process.calculateActions.stars.fluxType = "psfFlux_target"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.addSummaryPlot = False
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.yAxisLabel = "Output Mag - Ref Mag (mmag)"
        self.applyContext(VisitContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaPsfScatterPlot(TargetRefCatDeltaScatterPhotom):
    """Plot the difference in millimags between the PSF flux
    of a target catalog and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "{band}_psfFlux_target"


class TargetRefCatDeltaCModelScatterPlot(TargetRefCatDeltaScatterPhotom):
    """Plot the difference in millimags between the CModel flux
    of a target catalog and a reference catalog.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "{band}_cModelFlux_target"


class TargetRefCatDeltaPsfScatterVisitPlot(TargetRefCatDeltaScatterPhotomVisit):
    """Plot the difference in millimags between the PSF flux
    of a target catalog and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "psfFlux_target"


class TargetRefCatDeltaAp09ScatterVisitPlot(TargetRefCatDeltaScatterPhotomVisit):
    """Plot the difference in millimags between the aper 09 flux
    of a target catalog and a reference catalog.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "ap09Flux_target"


class TargetRefCatDeltaRAScatterPlot(TargetRefCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_target"
        )
        self.process.buildActions.yStars.buildAction.actionB = RAcosDec(raKey="ra_ref", decKey="dec_ref")

        self.produce.yAxisLabel = "RA$_{{target}}$ - RA$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaRAScatterVisitPlot(TargetRefCatDeltaScatterAstromVisit):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_target"
        )
        self.process.buildActions.yStars.buildAction.actionB = RAcosDec(raKey="ra_ref", decKey="dec_ref")

        self.produce.yAxisLabel = "RA$_{{target}}$ - RA$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaDecScatterVisitPlot(TargetRefCatDeltaScatterAstromVisit):
    """Plot the difference in milliseconds between the Decs of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.yStars.buildAction.actionB = LoadVector(vectorKey="dec_ref")

        self.produce.yAxisLabel = "RA$_{{target}}$ - RA$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaDecScatterPlot(TargetRefCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.yStars.buildAction.actionB = LoadVector(vectorKey="dec_ref")

        self.produce.yAxisLabel = "Dec$_{{target}}$ - Dec$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaSkyPlot(TargetRefCatDelta):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference of the quantity defined in the vector key parameter between
    the target and reference catalog as the color.
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


class TargetRefCatDeltaSkyPlotAstrom(TargetRefCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )

        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)

        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux_target"


class TargetRefCatDeltaSkyPlotAstromVisit(TargetRefCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars at
    the visit level, with the difference between the RA or Dec of
    the target and reference catalog as the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )

        self.applyContext(VisitContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaSkyPlotPhotomVisit(TargetRefCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the photometry of the target and reference catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = MagDiff()
        self.process.buildActions.zStars.col2 = "mag_ref"
        self.process.buildActions.zStars.fluxUnits2 = "mag(AB)"

        self.produce.plotName = "photomDiffSky"
        self.produce.zAxisLabel = "Output Mag - Ref Mag (mmag)"
        self.applyContext(VisitContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaSkyPlotPhotom(TargetRefCatDeltaSkyPlot):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the photometry of the target and reference catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars = MagDiff()
        self.process.buildActions.zStars.col2 = "{band}_mag_ref"
        self.process.buildActions.zStars.fluxUnits2 = "mag(AB)"

        self.produce.plotName = "photomDiffSky_{band}"
        self.produce.zAxisLabel = "Output Mag - Ref Mag (mmag)"
        self.applyContext(CoaddContext)
        self.applyContext(RefMatchContext)


class TargetRefCatDeltaPsfSkyPlot(TargetRefCatDeltaSkyPlotPhotom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the PSF photometry of the target and reference
    catalog as the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars.col1 = "{band}_psfFlux_target"


class TargetRefCatDeltaPsfSkyVisitPlot(TargetRefCatDeltaSkyPlotPhotomVisit):
    """Plot the RA/Dec distribution of stars, with the
    difference between the PSF photometry of the target and reference
    catalog as the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars.col1 = "psfFlux_target"


class TargetRefCatDeltaAp09SkyVisitPlot(TargetRefCatDeltaSkyPlotPhotomVisit):
    """Plot the RA/Dec distribution of stars, with the
    difference between the Ap09 photometry of the target and reference
    catalog as the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars.col1 = "ap09Flux_target"


class TargetRefCatDeltaCModelSkyPlot(TargetRefCatDeltaSkyPlotPhotom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the CModel photometry of the target and reference
    catalog as the color.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.zStars.col1 = "{band}_cModelFlux_target"


class TargetRefCatDeltaRASkyPlot(TargetRefCatDeltaSkyPlotAstrom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the RA of the target and reference catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.zStars.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_target"
        )
        self.process.buildActions.zStars.buildAction.actionB = RAcosDec(raKey="ra_ref", decKey="dec_ref")

        self.produce.plotName = "astromDiffSky_RA"
        self.produce.zAxisLabel = "RA$_{{target}}$ - RA$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaRASkyVisitPlot(TargetRefCatDeltaSkyPlotAstromVisit):
    """Plot the RA/Dec distribution of stars, with the
    difference between the RA of the target and reference catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.zStars.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_target"
        )
        self.process.buildActions.zStars.buildAction.actionB = RAcosDec(raKey="ra_ref", decKey="dec_ref")
        self.produce.plotName = "astromDiffSky_RA"
        self.produce.zAxisLabel = "RA$_{{target}}$ - RA$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaDecSkyVisitPlot(TargetRefCatDeltaSkyPlotAstromVisit):
    """Plot the RA/Dec distribution of stars, with the
    difference between the RA of the target and reference catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.zStars.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.zStars.buildAction.actionB = LoadVector(vectorKey="dec_ref")

        self.produce.plotName = "astromDiffSky_Dec"
        self.produce.zAxisLabel = "Dec$_{{target}}$ - Dec$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaDecSkyPlot(TargetRefCatDeltaSkyPlotAstrom):
    """Plot the RA/Dec distribution of stars, with the
    difference between the Dec of the target and reference catalog as
    the color.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.zStars.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.zStars.buildAction.actionB = LoadVector(vectorKey="dec_ref")

        self.produce.plotName = "astromDiffSky_Dec"
        self.produce.zAxisLabel = "Dec$_{{target}}$ - Dec$_{{ref}}$ (marcsec)"


class TargetRefCatDeltaMetrics(AnalysisTool):
    """Calculate the AA1 metric and the sigma MAD from the difference between
    the target and reference catalog coordinates.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        """Apply coadd options for the metrics. Applies the coadd plot flag
        selector and sets flux types.
        """
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness_target"

        self.applyContext(RefMatchContext)

        self.process.buildActions.mags.vectorKey = "{band}_psfFlux_target"

        self.produce.metric.newNames = {
            "AA1_RA": "{band}_AA1_RA_coadd",
            "AA1_sigmaMad_RA": "{band}_AA1_sigmaMad_RA_coadd",
            "AA1_Dec": "{band}_AA1_Dec_coadd",
            "AA1_sigmaMad_Dec": "{band}_AA1_sigmaMad_Dec_coadd",
            "AA1_tot": "{band}_AA1_tot_coadd",
            "AA1_sigmaMad_tot": "{band}_AA1_sigmaMad_tot_coadd",
        }

    def visitContext(self) -> None:
        """Apply visit options for the metrics. Applies the visit plot flag
        selector and sets flux types.
        """
        self.parameterizedBand = False
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()
        self.prep.selectors.starSelector.vectorKey = "extendedness_target"

        self.applyContext(RefMatchContext)

        self.process.buildActions.mags.vectorKey = "psfFlux_target"

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
        self.process.calculateActions.AA1_RA = MedianAction(vectorKey="brightStarsRA")
        self.process.calculateActions.AA1_sigmaMad_RA = SigmaMadAction(vectorKey="brightStarsRA")

        self.process.calculateActions.AA1_Dec = MedianAction(vectorKey="brightStarsDec")
        self.process.calculateActions.AA1_sigmaMad_Dec = SigmaMadAction(vectorKey="brightStarsDec")

        self.process.calculateActions.AA1_tot = MedianAction(vectorKey="brightStarsTot")
        self.process.calculateActions.AA1_sigmaMad_tot = SigmaMadAction(vectorKey="brightStarsTot")

        self.produce.metric.units = {
            "AA1_RA": "mas",
            "AA1_sigmaMad_RA": "mas",
            "AA1_Dec": "mas",
            "AA1_sigmaMad_Dec": "mas",
            "AA1_tot": "mas",
            "AA1_sigmaMad_tot": "mas",
        }
