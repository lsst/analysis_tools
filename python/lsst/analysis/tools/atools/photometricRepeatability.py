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
    "StellarPhotometricRepeatability",
    "StellarPhotometricResidualsFocalPlane",
)

from lsst.pex.config import Field

from ..actions.plot import FocalPlanePlot, HistPanel, HistPlot, HistStatsPanel
from ..actions.scalar.scalarActions import CountAction, FracThreshold, MedianAction, StdevAction
from ..actions.vector import (
    BandSelector,
    CalcSn,
    ConvertFluxToMag,
    LoadVector,
    MultiCriteriaDownselectVector,
    PerGroupStatistic,
    RangeSelector,
    ResidualWithPerGroupStatistic,
    SnSelector,
    ThresholdSelector,
)
from ..interfaces import AnalysisTool


class StellarPhotometricRepeatability(AnalysisTool):
    """Compute photometric repeatability from multiple measurements of a set of
    stars. First, a set of per-source quality criteria are applied. Second,
    the individual source measurements are grouped together by object index
    and per-group quantities are computed (e.g., a representative S/N for the
    group based on the median of associated per-source measurements). Third,
    additional per-group criteria are applied. Fourth, summary statistics are
    computed for the filtered groups.
    """

    fluxType = Field[str](doc="Flux type to calculate repeatability with", default="psfFlux")
    PA2Value = Field[float](
        doc="Used to compute the percent of individual measurements that deviate by more than PA2Value"
        "from the mean of each measurement (PF1). Units of PA2Value are mmag.",
        default=15.0,
    )

    def setDefaults(self):
        super().setDefaults()

        # Apply per-source selection criteria
        self.prep.selectors.bandSelector = BandSelector()

        # Compute per-group quantities
        self.process.buildActions.perGroupSn = PerGroupStatistic()
        self.process.buildActions.perGroupSn.buildAction = CalcSn()
        self.process.buildActions.perGroupSn.func = "median"
        self.process.buildActions.perGroupExtendedness = PerGroupStatistic()
        self.process.buildActions.perGroupExtendedness.buildAction.vectorKey = "extendedness"
        self.process.buildActions.perGroupExtendedness.func = "median"
        self.process.buildActions.perGroupCount = PerGroupStatistic()
        self.process.buildActions.perGroupCount.buildAction.vectorKey = f"{self.fluxType}"
        self.process.buildActions.perGroupCount.func = "count"
        # Use mmag units
        self.process.buildActions.perGroupStdev = PerGroupStatistic()
        self.process.buildActions.perGroupStdev.buildAction = ConvertFluxToMag(
            vectorKey=f"{self.fluxType}",
            returnMillimags=True,
        )
        self.process.buildActions.perGroupStdev.func = "std"

        # Filter on per-group quantities
        self.process.filterActions.perGroupStdevFiltered = MultiCriteriaDownselectVector(
            vectorKey="perGroupStdev"
        )
        self.process.filterActions.perGroupStdevFiltered.selectors.count = ThresholdSelector(
            vectorKey="perGroupCount",
            op="ge",
            threshold=3,
        )
        self.process.filterActions.perGroupStdevFiltered.selectors.sn = RangeSelector(
            vectorKey="perGroupSn",
            minimum=200,
        )
        self.process.filterActions.perGroupStdevFiltered.selectors.extendedness = ThresholdSelector(
            vectorKey="perGroupExtendedness",
            op="le",
            threshold=0.5,
        )

        # Compute summary statistics on filtered groups
        self.process.calculateActions.photRepeatStdev = MedianAction(vectorKey="perGroupStdevFiltered")
        self.process.calculateActions.photRepeatNsources = CountAction(vectorKey="perGroupStdevFiltered")

        self.produce.plot = HistPlot()

        self.produce.plot.panels["panel_rms"] = HistPanel()

        self.produce.plot.panels["panel_rms"].statsPanel = HistStatsPanel()
        self.produce.plot.panels["panel_rms"].statsPanel.statsLabels = ["N", "PA1", "PF1 %"]
        self.produce.plot.panels["panel_rms"].statsPanel.stat1 = ["photRepeatNsources"]
        self.produce.plot.panels["panel_rms"].statsPanel.stat2 = ["photRepeatStdev"]
        self.produce.plot.panels["panel_rms"].statsPanel.stat3 = ["photRepeatOutlier"]

        self.produce.plot.panels["panel_rms"].refRelativeToMedian = True

        self.produce.plot.panels["panel_rms"].label = "rms (mmag)"
        self.produce.plot.panels["panel_rms"].hists = dict(perGroupStdevFiltered="Filtered per group rms")

        self.produce.metric.units = {  # type: ignore
            "photRepeatStdev": "mmag",
            "photRepeatOutlier": "percent",
            "photRepeatNsources": "ct",
        }

    def finalize(self):
        super().finalize()
        self.process.buildActions.perGroupSn.buildAction.fluxType = f"{self.fluxType}"
        self.process.buildActions.perGroupCount.buildAction.vectorKey = f"{self.fluxType}"
        self.process.buildActions.perGroupStdev.buildAction = ConvertFluxToMag(
            vectorKey=f"{self.fluxType}",
            returnMillimags=True,
        )
        self.process.calculateActions.photRepeatOutlier = FracThreshold(
            vectorKey="perGroupStdevFiltered",
            op="ge",
            threshold=self.PA2Value,
            percent=True,
            relative_to_median=True,
        )

        if isinstance(self.produce.plot, HistPlot):
            self.produce.plot.panels["panel_rms"].referenceValue = self.PA2Value

        self.produce.metric.newNames = {
            "photRepeatStdev": "{band}_stellarPhotRepeatStdev",
            "photRepeatOutlier": "{band}_stellarPhotRepeatOutlierFraction",
            "photRepeatNsources": "{band}_ct",
        }


class StellarPhotometricResidualsFocalPlane(AnalysisTool):
    """Plot mean photometric residuals as a function of the position on the
    focal plane.

    First, a set of per-source quality criteria are applied. Second, the
    individual source measurements are grouped together by object index
    and the per-group magnitude is computed. The residuals between the
    individual sources and these magnitudes are then used to construct a plot
    showing the mean residual as a function of the focal-plane position.
    """

    fluxType = Field[str](doc="Flux type to calculate repeatability with", default="psfFlux")

    def setDefaults(self):
        super().setDefaults()

        # Apply per-source selection criteria
        self.prep.selectors.bandSelector = BandSelector()
        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 50

        self.process.buildActions.z = ResidualWithPerGroupStatistic()
        self.process.buildActions.z.buildAction = ConvertFluxToMag(
            vectorKey=f"{self.fluxType}",
            returnMillimags=True,
        )
        self.process.buildActions.z.func = "median"

        self.process.buildActions.x = LoadVector(vectorKey="x")
        self.process.buildActions.y = LoadVector(vectorKey="y")

        self.process.buildActions.detector = LoadVector(vectorKey="detector")

        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.threshold = 200
        self.process.buildActions.statMask.fluxType = "psfFlux"

        self.process.calculateActions.photResidFocalPlaneMedian = MedianAction(vectorKey="z")
        self.process.calculateActions.photResidFocalPlaneStdev = StdevAction(vectorKey="z")

        self.produce.plot = FocalPlanePlot()
        self.produce.plot.zAxisLabel = "Mag - Mag$_{mean}$ (mmag)"

        self.produce.metric.units = {  # type: ignore
            "photResidFocalPlaneStdev": "mmag",
            "photResidFocalPlaneMedian": "mmag",
        }
