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
    "CompletenessPurityTool",
    "TargetInjectedCatDelta",
    "TargetInjectedCatDeltaScatterAstrom",
    "TargetInjectedCatDeltaRAScatterPlot",
    "TargetInjectedCatDeltaDecScatterPlot",
    "TargetInjectedCatDeltaScatterPhotom",
    "TargetInjectedCatDeltaPsfScatterPlot",
    "TargetInjectedCatDeltaCModelScatterPlot",
    "TargetInjectedCatDeltaMetrics",
)

from lsst.pex.config import Field, ListField

from ..actions.plot.completenessPlot import CompletenessHist
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.scalar.scalarActions import MedianAction, SigmaMadAction
from ..actions.vector import (
    ConvertFluxToMag,
    ConvertUnits,
    DownselectVector,
    LoadVector,
    MagDiff,
    MagPercentileAction,
    RAcosDec,
    RangeSelector,
    SnSelector,
    SubtractVector,
)
from ..contexts import CoaddContext
from ..interfaces import AnalysisTool


class CompletenessPurityTool(AnalysisTool):
    """Plot the completeness or purity of injected sources by magnitude."""

    completenessPercentiles = ListField[float](doc="tmp", default=[16.0, 50.0, 84.0])

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.mag = ConvertFluxToMag(vectorKey="ref_{band}_flux")
        self.process.buildActions.matchDistance = LoadVector(vectorKey="match_distance")

        self.process.calculateActions.magPercentiles = MagPercentileAction(matchDistanceKey="match_distance", vectorKey="ref_{band}_flux", percentiles=self.completenessPercentiles)
        # self.process.calculateActions.magPercentiles.vectorKey = "ref_{band}_flux"
        # self.process.calculateActions.magPercentiles.matchDistanceKey = "match_distance"

        self.produce.plot = CompletenessHist()
        # self.produce.metric.units = {"mag20": "mag", "mag50": "mag"}
        # self.produce.metric.newNames = {"mag20": "{band}_mag20", "mag50": "{band}_mag50"}

    def finalize(self):
        for percentile in self.completenessPercentiles:
            name = f"mag{percentile}"
            # TODO: can't subscript. Try defining a dictionary and then appending or overwriting
            self.produce.metric.units[name] = "mag"
            self.produce.metric.newNames[name] = "{band}_" + name


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
        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 200
        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "{band}_psfFlux"
        self.process.calculateActions.stars.lowSNSelector.threshold = 300
        self.process.calculateActions.stars.highSNSelector.fluxType = "{band}_psfFlux"
        self.process.calculateActions.stars.highSNSelector.threshold = 2700
        self.process.calculateActions.stars.fluxType = "{band}_psfFlux"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"
        self.process.buildActions.patch = LoadVector()
        self.process.buildActions.patch.vectorKey = "patch"


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
        self.process.buildActions.xStars.vectorKey = "{band}_psfFlux"
        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.fluxType = "{band}_psfFlux"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.applyContext(CoaddContext)


class TargetInjectedCatDeltaRAScatterPlot(TargetInjectedCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the RA of a target catalog
    and an injected catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = RAcosDec(raKey="coord_ra", decKey="coord_dec")
        self.process.buildActions.yStars.buildAction.actionB = RAcosDec(raKey="ref_ra", decKey="ref_dec")

        self.produce.yAxisLabel = "RA$_{{output}}$ - RA$_{{input}}$ (marcsec)"


class TargetInjectedCatDeltaDecScatterPlot(TargetInjectedCatDeltaScatterAstrom):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and an injected catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.yStars.buildAction.actionA = LoadVector(vectorKey="coord_dec")
        self.process.buildActions.yStars.buildAction.actionB = LoadVector(vectorKey="ref_dec")

        self.produce.yAxisLabel = "Dec$_{{output}}$ - Dec$_{{input}}$ (marcsec)"


class TargetInjectedCatDeltaScatterPhotom(TargetInjectedCatDelta):
    """Plot the difference in millimags between a target catalog and an
    injected catalog for the flux type set in `setDefaults`.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars = MagDiff()
        self.process.buildActions.yStars.col2 = "ref_{band}_flux"

        self.process.buildActions.xStars = ConvertFluxToMag()
        self.process.buildActions.xStars.vectorKey = "{band}_psfFlux"

        self.produce = ScatterPlotWithTwoHists()
        self.produce.plotTypes = ["stars"]
        self.produce.magLabel = "PSF Magnitude (mag)"
        self.produce.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.yAxisLabel = "Output Mag - Input Mag (mmag)"
        self.applyContext(CoaddContext)


class TargetInjectedCatDeltaPsfScatterPlot(TargetInjectedCatDeltaScatterPhotom):
    """Plot the difference in millimags between the PSF flux
    of a target catalog and an injected catalog
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "{band}_psfFlux"


class TargetInjectedCatDeltaCModelScatterPlot(TargetInjectedCatDeltaScatterPhotom):
    """Plot the difference in millimags between the CModel flux
    of a target catalog and an injected catalog.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.yStars.col1 = "{band}_cModelFlux"


# TODO: get this to actually work for photometric metrics
class TargetInjectedCatDeltaMetrics(AnalysisTool):
    """Calculate the diff metric and the sigma MAD from the difference between
    the target and injected catalog coordinates and photometry.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        """Apply coadd options for the metrics. Applies the coadd plot flag
        selector and sets flux types.
        """

        # self.applyContext(RefMatchContext)

        self.process.buildActions.mags.vectorKey = "{band}_psfFlux"

        self.produce.metric.newNames = {
            "injected_RA_diff_median": "injected_{band}_RA_diff_median_coadd",
            "injected_RA_diff_sigmaMad": "injected_{band}_RA_diff_sigmaMad_coadd",
            "injected_Dec_diff_median": "injected_{band}_Dec_diff_median_coadd",
            "injected_Dec_diff_sigmaMad": "injected_{band}_Dec_diff_sigmaMad_coadd",
            "injected_photometry_diff_median": "injected_{band}_photometry_diff_median_coadd",
            "injected_photometry_diff_sigmaMad": "injected_{band}_photometry_diff_sigmaMad_coadd",
        }

    def setDefaults(self):
        super().setDefaults()

        # Calculate difference in RA
        self.process.buildActions.astromDiffRA = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )
        self.process.buildActions.astromDiffRA.buildAction.actionA = RAcosDec(
            raKey="coord_ra", decKey="coord_dec"
        )
        self.process.buildActions.astromDiffRA.buildAction.actionB = RAcosDec(
            raKey="ref_ra", decKey="ref_dec"
        )

        # Calculate difference in Dec
        self.process.buildActions.astromDiffDec = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )
        self.process.buildActions.astromDiffDec.buildAction.actionA = LoadVector(vectorKey="coord_dec")
        self.process.buildActions.astromDiffDec.buildAction.actionB = LoadVector(vectorKey="ref_dec")

        # Calculate difference in photometry
        self.process.buildActions.photDiff = MagDiff(col1="{band}_psfFlux", col2="ref_{band}_flux")

        # Filter down to only objects with mag 17-21.5
        self.process.buildActions.mags = ConvertFluxToMag()
        self.process.filterActions.brightStarsRA = DownselectVector(vectorKey="astromDiffRA")
        self.process.filterActions.brightStarsRA.selector = RangeSelector(
            vectorKey="mags", minimum=17, maximum=21.5
        )
        self.process.filterActions.brightStarsDec = DownselectVector(vectorKey="astromDiffDec")
        self.process.filterActions.brightStarsDec.selector = RangeSelector(
            vectorKey="mags", minimum=17, maximum=21.5
        )
        self.process.filterActions.brightStarsPhot = DownselectVector(vectorKey="photDiff")
        self.process.filterActions.brightStarsPhot.selector = RangeSelector(
            vectorKey="mags", minimum=17, maximum=21.5
        )

        # Calculate median and sigmaMad
        self.process.calculateActions.injected_RA_diff_median = MedianAction(vectorKey="brightStarsRA")
        self.process.calculateActions.injected_RA_diff_sigmaMad = SigmaMadAction(vectorKey="brightStarsRA")

        self.process.calculateActions.injected_Dec_diff_median = MedianAction(vectorKey="brightStarsDec")
        self.process.calculateActions.injected_Dec_diff_sigmaMad = SigmaMadAction(vectorKey="brightStarsDec")

        self.process.calculateActions.injected_phot_diff_median = MedianAction(vectorKey="brightStarsPhot")
        self.process.calculateActions.injected_phot_diff_sigmaMad = SigmaMadAction(
            vectorKey="brightStarsPhot"
        )

        self.produce.metric.units = {
            "injected_RA_diff_median": "mas",
            "injected_RA_diff_sigmaMad": "mas",
            "injected_Dec_diff_median": "mas",
            "injected_Dec_diff_sigmaMad": "mas",
            "injected_phot_diff_median": "mmag",
            "injected_phot_diff_sigmaMad": "mmag",
        }
