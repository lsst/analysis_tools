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
)

from lsst.pex.config import Field

from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import (
    ConvertFluxToMag,
    ConvertUnits,
    DownselectVector,
    LoadVector,
    RAcosDec,
    SnSelector,
    StarSelector,
    SubtractVector,
    VectorSelector,
)
from ..interfaces import AnalysisTool
from .genericPrep import CoaddPrep, VisitPrep


class TargetRefCatDelta(AnalysisTool):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        self.prep = CoaddPrep()
        self.process.buildActions.starSelector.vectorKey = "{band}_extendedness"
        self.process.buildActions.mags = ConvertFluxToMag(vectorKey="{band}_psfFlux")
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="{band}_psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="{band}_psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

    def visitContext(self) -> None:
        self.parameterizedBand = False
        self.prep = VisitPrep()
        self.process.buildActions.starSelector.vectorKey = "extendedness"
        self.process.buildActions.mags = ConvertFluxToMag(vectorKey="psfFlux")
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

    def setDefaults(self, coordinate):
        super().setDefaults()

        self.process.buildActions.starSelector = StarSelector()
        self.process.buildActions.astromDiff = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )

        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="astromDiff", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.stars = ScatterPlotStatsAction(vectorKey="yStars")
        self.process.calculateActions.stars.lowSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.highSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.fluxType = "psfFlux"

        self.produce.plot = ScatterPlotWithTwoHists()

        self.produce.plot.plotTypes = ["stars"]
        self.produce.plot.xAxisLabel = "PSF Magnitude (mag)"
        self.produce.plot.yAxisLabel = f"${coordinate}_{{target}} - {coordinate}_{{ref}}$ (marcsec)"
        self.produce.plot.magLabel = "PSF Magnitude (mag)"


class TargetRefCatDeltaRAScatterPlot(TargetRefCatDelta):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(coordinate="RA")
        self.process.buildActions.astromDiff.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_ref"
        )
        self.process.buildActions.astromDiff.buildAction.actionB = RAcosDec(
            raKey="coord_ra_ref", decKey="coord_dec_ref"
        )


class TargetRefCatDeltaDecScatterPlot(TargetRefCatDelta):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults(coordinate="Dec")
        self.process.buildActions.astromDiff.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")
        self.process.buildActions.astromDiff.buildAction.actionB = LoadVector(vectorKey="coord_dec_ref")


class TargetRefCatDeltaSkyPlot(AnalysisTool):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def coaddContext(self) -> None:
        self.prep = CoaddPrep()

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"

    def visitContext(self) -> None:
        self.parameterizedBand = False
        self.prep = VisitPrep()

        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "psfFlux"

    def setDefaults(self, coordinate):
        super().setDefaults()

        self.process.buildActions.zStars = ConvertUnits(
            buildAction=SubtractVector, inUnit="degree", outUnit="milliarcsecond"
        )
        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra_target"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec_target"
        self.process.buildActions.patch = LoadVector()
        self.process.buildActions.patch.vectorKey = "patch"

        self.produce = SkyPlot()
        self.produce.plotTypes = ["stars"]
        self.produce.plotName = f"astromDiffSky_{coordinate}"
        self.produce.xAxisLabel = "R.A. (degrees)"
        self.produce.yAxisLabel = "Dec. (degrees)"
        self.produce.zAxisLabel = f"${coordinate}_{{target}} - {coordinate}_{{ref}}$ (marcsec)"
        self.produce.plotOutlines = False


class TargetRefCatDeltaRASkyPlot(TargetRefCatDeltaSkyPlot):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults(coordinate="RA")
        self.process.buildActions.zStars.buildAction.actionA = RAcosDec(
            raKey="coord_ra_target", decKey="coord_dec_ref"
        )
        self.process.buildActions.zStars.buildAction.actionB = RAcosDec(
            raKey="coord_ra_ref", decKey="coord_dec_ref"
        )


class TargetRefCatDeltaDecSkyPlot(TargetRefCatDeltaSkyPlot):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults(coordinate="Dec")

        self.process.buildActions.zStars.buildAction.actionA = LoadVector(vectorKey="coord_dec_target")

        self.process.buildActions.zStars.buildAction.actionB = LoadVector(vectorKey="coord_dec_ref")
