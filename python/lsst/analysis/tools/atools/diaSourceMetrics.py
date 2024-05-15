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
    "NumDiaSourcesMetric",
    "NumDipolesMetric",
    "NumDiaSourcesSelectionMetric",
    "DiaSourcesGoodVsBadRatioMetric",
)

from lsst.pex.config import Field

from ..actions.scalar import CountAction, DivideScalar
from ..actions.vector import DownselectVector, FlagSelector, GoodDiaSourceSelector
from ..interfaces import AnalysisTool
from ..contexts import DrpContext


class NumDiaSourcesMetric(AnalysisTool):
    """Calculate the number of DIA Sources that do not have known
    bad/quality flags set to true.
    """

    def setDefaults(self):
        super().setDefaults()

        # select dia sources that do not have bad flags
        #self.prep.selectors.goodDiaSourceSelector = GoodDiaSourceSelector()

        # Count the number of dia sources left after filtering
        self.process.filterActions.goodDiaSources = DownselectVector(
            vectorKey="diaSourceId",
            selector=GoodDiaSourceSelector()
        )
        self.process.calculateActions.numGoodDiaSources = CountAction(vectorKey="goodDiaSources")
        self.process.calculateActions.numDiaSources = CountAction(vectorKey="diaSourceId")
        self.process.calculateActions.diaSourcesGoodVsAllRatio = DivideScalar(
            actionA=self.process.calculateActions.numGoodDiaSources,
            actionB=self.process.calculateActions.numDiaSources,
        )
        # the units for the quantity (count, an astropy quantity)
        #self.produce.metric.units = {"numDiaSources": "ct",
        #                             "diaSourcesGoodVsAllRatio": "ct"
        #                             }
        self.produce.metric.units = {"diaSourcesGoodVsAllRatio": "ct"}
        self.applyContext(DrpContext)


class NumDipolesMetric(AnalysisTool):
    """Calculate the number of dipoles with NaN values excluded."""

    def setDefaults(self):
        super().setDefaults()

        # select all diaSources flagged as dipole
        self.prep.selectors.flags = FlagSelector(selectWhenTrue=["isDipole"])

        # count the number of dipoles
        self.process.buildActions.numDipoles = CountAction(vectorKey="isDipole")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"numDipoles": "ct"}


class NumDiaSourcesSelectionMetric(AnalysisTool):
    """Count the number of DIA Sources for a given threshold."""

    metricName = Field[str](doc="Name to use for output metric")

    def setDefaults(self):
        super().setDefaults()

        # Count dia sources with reliability lower than the threshold
        self.process.calculateActions.countingAction = CountAction

        # The units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"countingAction": "ct"}

    def finalize(self):
        self.produce.metric.newNames = {"countingAction": self.metricName}


class DiaSourcesGoodVsBadRatioMetric(AnalysisTool):
    """Calculate the ratio of 'good' vs 'bad' DIA Sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count dia sources with reliability higher than the threshold
        self.process.buildActions.numDiaSourcesHighReliability = CountAction(
            op="gt", threshold=0.9, vectorKey="reliability"
        )

        # Count dia sources with reliability lower than the threshold
        self.process.buildActions.numDiaSourcesLowReliability = CountAction(
            op="lt", threshold=0.1, vectorKey="reliability"
        )

        # Calculate ratio of good vs bad DIA Sources
        self.process.calculateActions.DiaSourcesGoodVsBadRatio = DivideScalar(
            actionA=self.process.buildActions.numDiaSourcesHighReliability,
            actionB=self.process.buildActions.numDiaSourcesLowReliability,
        )

        # The units for the quantity (dimensionless, an astropy quantity)
        self.produce.metric.units = {"DiaSourcesGoodVsBadRatio": ""}
