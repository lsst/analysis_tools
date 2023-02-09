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

__all__ = ("SourcesTool",)

from ..actions.plot.barPlots import BarPanel, BarPlot
from ..actions.scalar import CountAction, CountUniqueAction, FracThreshold
from ..actions.vector import DownselectVector, LoadVector, ThresholdSelector, VectorSelector
from ..interfaces import AnalysisTool


class SourcesTool(AnalysisTool):
    """Plot a histogram of the associated and unassociated sources."""

    def setDefaults(self, **kwargs):
        super().setDefaults()

        self.process.buildActions.loadVectorSources = LoadVector()
        self.process.buildActions.associatedVectorSelector = ThresholdSelector()
        self.process.buildActions.unassociatedVectorSelector = ThresholdSelector()

        # assign keys for PSF and AP Flux
        self.process.buildActions.loadVectorSources.vectorKey = "nDiaSources"
        self.process.buildActions.associatedVectorSelector.vectorKey = "nDiaSources"
        self.process.buildActions.associatedVectorSelector.op = "gt"
        self.process.buildActions.associatedVectorSelector.threshold = 1.0
        self.process.buildActions.unassociatedVectorSelector.vectorKey = "nDiaSources"
        self.process.buildActions.unassociatedVectorSelector.op = "le"
        self.process.buildActions.unassociatedVectorSelector.threshold = 1.0

        self.process.filterActions.allSources = VectorSelector(vectorKey="nDiaSources")
        self.process.filterActions.associatedVector = DownselectVector(
            vectorKey="nDiaSources", selector=self.process.buildActions.associatedVectorSelector
        )
        self.process.filterActions.unassociatedVector = DownselectVector(
            vectorKey="nDiaSources", selector=self.process.buildActions.unassociatedVectorSelector
        )

        self.process.buildActions.uniqueSources = CountUniqueAction(vectorKey="nDiaSources")
        self.process.calculateActions.associatedCount = CountAction(vectorKey="associatedVector")
        self.process.calculateActions.unassociatedCount = CountAction(vectorKey="unassociatedVector")

        self.process.calculateActions.associatedPercent = FracThreshold(
            op="gt", threshold=1.0, vectorKey="nDiaSources"
        )
        self.process.calculateActions.unassociatedPercent = FracThreshold(
            op="le", threshold=1.0, vectorKey="nDiaSources"
        )

        self.process.calculateActions.associatedCount = CountAction(vectorKey="associatedVector")
        self.process.calculateActions.unassociatedCount = CountAction(vectorKey="unassociatedVector")

        self.produce.plot = BarPlot()
        self.produce.plot.panels["panel_source"] = BarPanel()
        self.produce.plot.panels["panel_source"].label = "N Assoc and Unassoc Sources"
        self.produce.plot.panels["panel_source"].bars = dict(
            associatedVector="Associated Sources",
            unassociatedVector="Unassociated Sources",
        )

        self.produce.metric.units = {
            "associatedPercent": "count",
            "unassociatedPercent": "count",
            "associatedCount": "count",
            "unassociatedCount": "count",
            "uniqueSources": "count",
        }
