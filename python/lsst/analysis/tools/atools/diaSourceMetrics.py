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
    "NumDiaSourcesAllMetric",
    "NumDiaSourcesMetric",
    "NumDipolesMetric",
    "NumStreakSourcesMetric"
)

from ..actions.scalar import CountAction
from ..actions.vector import FlagSelector, GoodDiaSourceSelector
from ..interfaces import AnalysisTool


class NumDiaSourcesAllMetric(AnalysisTool):
    """Calculate the number of DIA Sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.numDiaSourcesAll = CountAction(vectorKey="diaSourceId")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"numDiaSourcesAll": "ct"}


class NumDiaSourcesMetric(AnalysisTool):
    """Calculate the number of DIA Sources that do not have known
    bad/quality flags set to true.
    """

    def setDefaults(self):
        super().setDefaults()

        # select dia sources that do not have bad flags
        self.prep.selectors.goodDiaSourceSelector = GoodDiaSourceSelector()

        # Count the number of dia sources left after filtering
        self.process.calculateActions.numDiaSources = CountAction(vectorKey="diaSourceId")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"numDiaSources": "ct"}


class NumDipolesMetric(AnalysisTool):
    """Calculate the number of dipoles."""

    def setDefaults(self):
        super().setDefaults()

        # select all diaSources flagged as dipole
        self.prep.selectors.flags = FlagSelector(selectWhenTrue=["isDipole"])

        # count the number of dipoles
        self.process.buildActions.numDipoles = CountAction(vectorKey="isDipole")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"numDipoles": "ct"}


class NumStreakSourcesMetric(AnalysisTool):
    """Calculate the number of diaSources with a streak in their footprint."""

    def setDefaults(self):
        super().setDefaults()

        # Select all diaSources with the streak mask set
        self.prep.selectors.flags = FlagSelector(selectWhenTrue=["base_PixelFlags_flag_streak"])

        # Count the number of diaSources in a streak
        self.process.buildActions.numStreakSources = CountAction(vectorKey="base_PixelFlags_flag_streak")

        # Units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"numStreakSources": "ct"}
