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
    "AstrometricRelativeRepeatability",
    "StellarAstrometricResidualsRAFocalPlanePlot",
    "StellarAstrometricResidualsDecFocalPlanePlot",
    "StellarAstrometricResidualStdDevRAFocalPlanePlot",
    "StellarAstrometricResidualStdDevDecFocalPlanePlot",
    "StellarAstrometricResidualsRASkyPlot",
    "StellarAstrometricResidualsDecSkyPlot",
)

from lsst.pex.config import Field

from ..actions.keyedData import CalcRelativeDistances
from ..actions.plot import FocalPlanePlot, HistPanel, HistPlot, SkyPlot
from ..actions.vector import (
    BandSelector,
    ConvertFluxToMag,
    ConvertUnits,
    DownselectVector,
    LoadVector,
    RAcosDec,
    RangeSelector,
    ResidualWithPerGroupStatistic,
    SnSelector,
    ThresholdSelector,
)
from ..interfaces import AnalysisTool


class StellarAstrometricResidualsBase(AnalysisTool):
    """Plot mean astrometric residuals.

    The individual source measurements are grouped together by object index
    and the per-group centroid is computed. The residuals between the
    individual sources and these centroids are then used to construct a plot
    showing the mean residual as a function of the focal-plane or sky position.
    """

    def setDefaults(self):
        super().setDefaults()

        # Apply per-source selection criteria
        self.prep.selectors.bandSelector = BandSelector()

        self.process.buildActions.mags = ConvertFluxToMag(vectorKey="psfFlux")
        self.process.buildActions.residual = ConvertUnits()
        self.process.buildActions.residual.inUnit = "degree"
        self.process.buildActions.residual.outUnit = "marcsec"

        self.process.buildActions.residual.buildAction = ResidualWithPerGroupStatistic()

        self.process.buildActions.x = LoadVector(vectorKey="x")
        self.process.buildActions.y = LoadVector(vectorKey="y")

        self.process.buildActions.detector = LoadVector(vectorKey="detector")

        self.process.filterActions.x = DownselectVector(vectorKey="x")
        self.process.filterActions.x.selector = ThresholdSelector(
            vectorKey="mags",
            op="le",
            threshold=24,
        )
        self.process.filterActions.y = DownselectVector(
            vectorKey="y", selector=self.process.filterActions.x.selector
        )
        self.process.filterActions.z = DownselectVector(
            vectorKey="residual", selector=self.process.filterActions.x.selector
        )
        self.process.filterActions.detector = DownselectVector(
            vectorKey="detector", selector=self.process.filterActions.x.selector
        )

        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.threshold = 0
        self.process.buildActions.statMask.fluxType = "psfFlux"


class StellarAstrometricResidualsRASkyPlot(StellarAstrometricResidualsBase):
    """Plot mean astrometric residuals in RA as a function of the position in
    RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults()

        # Compute per-group quantities
        self.process.buildActions.residual.buildAction.buildAction = RAcosDec()
        self.process.buildActions.x = LoadVector(vectorKey="coord_ra")
        self.process.buildActions.y = LoadVector(vectorKey="coord_dec")

        self.produce = SkyPlot()

        self.produce.plotTypes = ["any"]
        self.produce.plotName = "ra_residuals"
        self.produce.xAxisLabel = "R.A. (degrees)"
        self.produce.yAxisLabel = "Dec. (degrees)"
        self.produce.zAxisLabel = "RAcos(Dec) - RAcos(Dec)$_{mean}$"


class StellarAstrometricResidualsDecSkyPlot(StellarAstrometricResidualsBase):
    """Plot mean astrometric residuals in RA as a function of the position in
    RA and Dec.
    """

    def setDefaults(self):
        super().setDefaults()

        # Compute per-group quantities
        self.process.buildActions.residual.buildAction.buildAction.vectorKey = "coord_dec"
        self.process.buildActions.x = LoadVector(vectorKey="coord_ra")
        self.process.buildActions.y = LoadVector(vectorKey="coord_dec")

        self.produce = SkyPlot()

        self.produce.plotTypes = ["any"]
        self.produce.plotName = "ra_residuals"
        self.produce.xAxisLabel = "R.A. (degrees)"
        self.produce.yAxisLabel = "Dec. (degrees)"
        self.produce.zAxisLabel = "RAcos(Dec) - RAcos(Dec)$_{mean}$"


class StellarAstrometricResidualsRAFocalPlanePlot(StellarAstrometricResidualsBase):
    """Plot mean astrometric residuals in RA as a function of the focal plane
    position.
    """

    def setDefaults(self):
        super().setDefaults()

        # Compute per-group quantities
        self.process.buildActions.residual.buildAction.buildAction = RAcosDec()

        self.produce = FocalPlanePlot()
        self.produce.zAxisLabel = "RAcos(Dec) - RAcos(Dec)$_{mean}$ (mArcsec)"


class StellarAstrometricResidualStdDevRAFocalPlanePlot(StellarAstrometricResidualsBase):
    """Plot mean astrometric residuals in RA as a function of the focal plane
    position.
    """

    def setDefaults(self):
        super().setDefaults()

        # Compute per-group quantities
        self.process.buildActions.residual.buildAction.buildAction = RAcosDec()

        self.produce = FocalPlanePlot()
        self.produce.statistic = "std"
        self.produce.zAxisLabel = "Std(RAcos(Dec) - RAcos(Dec)$_{mean}$) (mArcsec)"


class StellarAstrometricResidualsDecFocalPlanePlot(StellarAstrometricResidualsBase):
    """Plot mean astrometric residuals in RA as a function of the focal plane
    position.
    """

    def setDefaults(self):
        super().setDefaults()

        # Compute per-group quantities
        self.process.buildActions.residual.buildAction.buildAction.vectorKey = "coord_dec"

        self.produce = FocalPlanePlot()
        self.produce.zAxisLabel = "Dec - Dec$_{mean}$ (mArcsec)"


class StellarAstrometricResidualStdDevDecFocalPlanePlot(StellarAstrometricResidualsBase):
    """Plot mean astrometric residuals in RA as a function of the focal plane
    position.
    """

    def setDefaults(self):
        super().setDefaults()

        # Compute per-group quantities
        self.process.buildActions.residual.buildAction.buildAction.vectorKey = "coord_dec"

        self.produce = FocalPlanePlot()
        self.produce.statistic = "std"
        self.produce.zAxisLabel = "Std(Dec - Dec$_{mean}$) (mArcsec)"


class AstrometricRelativeRepeatability(AnalysisTool):
    """Calculate the AMx, ADx, AFx metrics and make histograms showing the data
    used to compute the metrics.
    """

    fluxType = Field[str](doc="Flux type to calculate repeatability with", default="psfFlux")
    xValue = Field[int](doc="Metric suffix corresponding to annulus size (1, 2, or 3)", default=1)

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.bandSelector = BandSelector()
        # Following what was done in faro, only sources with S/N between 50
        # and 50000 are included. The other filtering that was done in faro
        # is now covered by only including sources from isolated_star_sources.
        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.threshold = 50
        self.prep.selectors.snSelector.maxSN = 50000

        # Select only sources with magnitude between 17 and 21.5
        self.process.filterActions.coord_ra = DownselectVector(vectorKey="coord_ra")
        self.process.filterActions.coord_ra.selector = RangeSelector(
            vectorKey="mags", minimum=17, maximum=21.5
        )
        self.process.filterActions.coord_dec = DownselectVector(
            vectorKey="coord_dec", selector=self.process.filterActions.coord_ra.selector
        )
        self.process.filterActions.obj_index = DownselectVector(
            vectorKey="obj_index", selector=self.process.filterActions.coord_ra.selector
        )
        self.process.filterActions.visit = DownselectVector(
            vectorKey="visit", selector=self.process.filterActions.coord_ra.selector
        )

        self.process.calculateActions.rms = CalcRelativeDistances()

        self.produce.metric.units = {
            "AMx": "mas",
            "AFx": "percent",
            "ADx": "mas",
        }

        self.produce.plot = HistPlot()

        self.produce.plot.panels["panel_sep"] = HistPanel()
        self.produce.plot.panels["panel_sep"].hists = dict(separationResiduals="Source separations")
        self.produce.plot.panels["panel_sep"].label = "Separation Distances (marcsec)"

        self.produce.plot.panels["panel_rms"] = HistPanel()
        self.produce.plot.panels["panel_rms"].hists = dict(rmsDistances="Object RMS")
        self.produce.plot.panels["panel_rms"].label = "Per-Object RMS (marcsec)"
        # TODO: DM-39163 add reference lines for ADx, AMx, and AFx.

    def finalize(self):
        super().finalize()
        self.prep.selectors.snSelector.fluxType = self.fluxType
        self.process.buildActions.mags = ConvertFluxToMag(vectorKey=self.fluxType)

        self.produce.metric.newNames = {
            "AMx": f"{{band}}_AM{self.xValue}",
            "AFx": f"{{band}}_AF{self.xValue}",
            "ADx": f"{{band}}_AD{self.xValue}",
        }
