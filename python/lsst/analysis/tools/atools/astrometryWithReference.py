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
    "CoordinateConfig",
    "TargetRefCatDeltaRAScatterPlot",
    "TargetRefCatDeltaDecScatterPlot",
    "TargetRefCatDeltaRASkyPlot",
    "TargetRefCatDeltaDecSkyPlot",
)

from lsst.pex.config import ChoiceField, Config, Field

from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import (
    AstromDiff,
    DownselectVector,
    LoadVector,
    MagColumnNanoJansky,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from ..interfaces import AnalysisTool
from .coaddVisit import CoaddVisitConfig
from .genericPrep import CoaddPrep, VisitPrep


class CoordinateConfig(Config):
    coordinate = ChoiceField[str](
        doc="The name of the sky coordinate",
        allowed={"RA": "Right Ascension", "Dec": "Declination"},
        optional=False,
    )


class TargetRefCatDelta(AnalysisTool, CoaddVisitConfig, CoordinateConfig):
    """Plot the difference in milliseconds between a target catalog and a
    reference catalog for the coordinate set in `setDefaults`.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def finalize(self):
        AnalysisTool().finalize()
        if self.context == "coadd":
            self.prep = CoaddPrep()
            self.process.buildActions.starSelector.vectorKey = "{band}_extendedness"
            self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="{band}_psfFlux")
            self.process.filterActions.psfFlux = DownselectVector(
                vectorKey="{band}_psfFlux", selector=VectorSelector(vectorKey="starSelector")
            )
            self.process.filterActions.psfFluxErr = DownselectVector(
                vectorKey="{band}_psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
            )
        elif self.context == "visit":
            self.parameterizedBand = False
            self.prep = VisitPrep()
            self.process.buildActions.starSelector.vectorKey = "extendedness"
            self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="psfFlux")
            self.process.filterActions.psfFlux = DownselectVector(
                vectorKey="psfFlux", selector=VectorSelector(vectorKey="starSelector")
            )
            self.process.filterActions.psfFluxErr = DownselectVector(
                vectorKey="psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
            )
        else:
            raise ValueError(f"Unsupported {self.context=}")

        coordStr = self.coordinate.lower()
        self.process.buildActions.astromDiff = AstromDiff(
            col1=f"coord_{coordStr}_target", col2=f"coord_{coordStr}_ref"
        )
        self.produce.plot.yAxisLabel = f"${self.coordinate}_{{target}} - {self.coordinate}_{{ref}}$ (marcsec)"

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.starSelector = StarSelector()

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
        # Placeholder
        self.produce.plot.yAxisLabel = ""
        self.produce.plot.magLabel = "PSF Magnitude (mag)"


class TargetRefCatDeltaRAScatterPlot(TargetRefCatDelta):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.coordinate = "RA"


class TargetRefCatDeltaDecScatterPlot(TargetRefCatDelta):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog
    """

    def setDefaults(self):
        super().setDefaults()
        self.coordinate = "Dec"


class TargetRefCatDeltaSkyPlot(AnalysisTool, CoaddVisitConfig, CoordinateConfig):
    """Base class for plotting the RA/Dec distribution of stars, with the
    difference between the RA or Dec of the target and reference catalog as
    the color.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=True
    )

    def finalize(self):
        AnalysisTool().finalize()
        if self.context == "coadd":
            self.prep = CoaddPrep()

            self.process.buildActions.starStatMask = SnSelector()
            self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"
        elif self.context == "visit":
            self.parameterizedBand = False
            self.prep = VisitPrep()

            self.process.buildActions.starStatMask = SnSelector()
            self.process.buildActions.starStatMask.fluxType = "psfFlux"
        else:
            raise ValueError(f"Unsupported {self.context=}")

        coordStr = self.coordinate.lower()
        self.process.buildActions.zStars = AstromDiff(
            col1=f"coord_{coordStr}_target", col2=f"coord_{coordStr}_ref"
        )
        self.produce.plotName = f"astromDiffSky_{self.coordinate}"
        self.produce.zAxisLabel = f"${self.coordinate}_{{target}} - {self.coordinate}_{{ref}}$ (marcsec)"

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
        # Placeholders
        self.produce.zAxisLabel = ""
        self.produce.plotName = ""
        self.produce.plotOutlines = False


class TargetRefCatDeltaRASkyPlot(TargetRefCatDeltaSkyPlot):
    """Plot the difference in milliseconds between the RA of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults()
        self.coordinate = "RA"


class TargetRefCatDeltaDecSkyPlot(TargetRefCatDeltaSkyPlot):
    """Plot the difference in milliseconds between the Dec of a target catalog
    and a reference catalog as a function of RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults()
        self.coordinate = "Dec"
