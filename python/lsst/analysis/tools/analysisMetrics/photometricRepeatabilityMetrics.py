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

__all__ = ("StellarPhotometricRepeatabilityMetric",)

from ..actions.scalar.scalarActions import FracThreshold, MedianAction
from ..actions.vector import (
    BandSelector,
    DownselectVector,
    MagColumnNanoJansky,
    PerGroupStatistic,
    Sn,
    ThresholdSelector,
)
from ..interfaces import AnalysisMetric


class StellarPhotometricRepeatabilityMetric(AnalysisMetric):
    """Compute photometric repeatability from multiple observations of the same set of objects."""

    fluxType: str = "psfFlux"

    def setDefaults(self):
        super().setDefaults()

        # Apply per-source selection criteria
        self.prep.selectors.bandSelector = BandSelector()

        # Compute per-group quantities
        self.process.buildActions.perGroupSn = PerGroupStatistic()
        self.process.buildActions.perGroupSn.buildAction = Sn(fluxType=f"{self.fluxType}")
        self.process.buildActions.perGroupSn.func = "median"
        self.process.buildActions.perGroupExtendedness = PerGroupStatistic()
        self.process.buildActions.perGroupExtendedness.buildAction.vectorKey = "extendedness"
        self.process.buildActions.perGroupExtendedness.func = "median"
        self.process.buildActions.perGroupCount = PerGroupStatistic()
        self.process.buildActions.perGroupCount.buildAction.vectorKey = f"{self.fluxType}"
        self.process.buildActions.perGroupStdev = PerGroupStatistic()
        self.process.buildActions.perGroupStdev.buildAction = MagColumnNanoJansky(
            vectorKey=f"{self.fluxType}"
        )
        self.process.buildActions.perGroupStdev.func = "std"

        # Filter on per-group quantities
        """
        self.process.filterActions.perGroupStdevFiltered = DownselectVector(vectorKey="perGroupStdev")
        self.process.filterActions.perGroupStdevFiltered.selector = AndSelector()
        self.process.filterActions.perGroupStdevFiltered.selector.selectors.count = ThresholdSelector(
            vectorKey="perGroupCount", op="ge", threshold=3,
        )
        self.process.filterActions.perGroupStdevFiltered.selector.selectors.sn = ThresholdSelector(
            vectorKey="perGroupSn", op="ge", threshold=200,
        )
        self.process.filterActions.perGroupStdevFiltered.selector.selectors.extendedness = ThresholdSelector(
            vectorKey="perGroupExtendedness", op="le", threshold=0.5,
        )
        """
        self.process.filterActions.perGroupStdevFiltered = DownselectVector(vectorKey="perGroupStdev")
        self.process.filterActions.perGroupStdevFiltered.selectors.count = ThresholdSelector(
            vectorKey="perGroupCount",
            op="ge",
            threshold=3,
        )
        self.process.filterActions.perGroupStdevFiltered.selectors.sn = ThresholdSelector(
            vectorKey="perGroupSn",
            op="ge",
            threshold=200,
        )
        self.process.filterActions.perGroupStdevFiltered.selectors.extendedness = ThresholdSelector(
            vectorKey="perGroupExtendedness",
            op="le",
            threshold=0.5,
        )
        """
        self.process.filterActions.perGroupStdevFiltered = KeyedDataSelectorAction(vectorKeys=["perGroupStdev"])
        self.process.filterActions.perGroupStdevFiltered.selectors.count = ThresholdSelector(
            vectorKey="perGroupCount", op="ge", threshold=3,
        )
        self.process.filterActions.perGroupStdevFiltered.selectors.sn = ThresholdSelector(
            vectorKey="perGroupSn", op="ge", threshold=200,
        )
        self.process.filterActions.perGroupStdevFiltered.selectors.extendedness = ThresholdSelector(
            vectorKey="perGroupExtendedness", op="le", threshold=0.5,
        )
        """

        self.process.calculateActions.photRepeatStdev = MedianAction(vectorKey="perGroupStdevFiltered")
        self.process.calculateActions.photRepeatOutlier = FracThreshold(
            vectorKey="perGroupStdevFiltered",
            op="ge",
            threshold=0.015,
            percent=True,
        )

        self.produce.units = {  # type: ignore
            "photRepeatStdev": "mag",
            "photRepeatOutlier": "percent",
        }
        self.produce.newNames = {
            "photRepeatStdev": "{band}_stellarPhotRepeatStdev",
            "photRepeatOutlier": "{band}_stellarPhotRepeatOutlierFraction",
        }
