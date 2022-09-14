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

__all__ = ("NumSsObjectsMetric",)

from ..actions.scalar import CountAction
from ..actions.vector import DownselectVector, ThresholdSelector
from ..interfaces import AnalysisMetric


class NumSsObjectsMetric(AnalysisMetric):
    """Count the number of DIASources that are associated with known
    Solar System Objects.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.thresholdSelector = ThresholdSelector()

        # select nonzero ssObjectId
        self.process.buildActions.thresholdSelector.vectorKey = "ssObjectId"
        self.process.buildActions.thresholdSelector.op = "ge"
        self.process.buildActions.thresholdSelector.threshold = 1

        # the final name in the qualification is used as a key to insert
        # the calculation into KeyedData
        self.process.filterActions.allSSOs = DownselectVector(
            vectorKey="ssObjectId", selector=self.process.buildActions.thresholdSelector
        )

        self.process.calculateActions.NumSsObjectsMetric = CountAction(vectorKey="allSSOs")

        self.produce.units = {"NumSsObjectsMetric": "ct"}
