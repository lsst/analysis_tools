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

__all__ = (
    "NumDiaSourcesMetric",
    "NumStreakDiaSourcesMetric",
    "NumStreakCenterDiaSourcesMetric",
)

from ..actions.scalar import CountAction
from ..actions.vector import FlagSelector, GoodDiaSourceSelector
from ..interfaces import AnalysisTool
from ..contexts import DrpContext


class NumDiaSourcesMetric(AnalysisTool):
    """Count all DiaSources that do not have known bad/quality flags."""

    def setDefaults(self):
        super().setDefaults()

        # Select "good" DiaSources that are not obvious garbage
        self.prep.selectors.goodDiaSourceSelector = GoodDiaSourceSelector()
        # Count them
        self.process.calculateActions.numDiaSources = CountAction(vectorKey="parentDiaSourceId")
        # Set the units for the resulting astropy quantity (ct is count)
        self.produce.metric.units = {"numDiaSources": "ct"}
        # Use, e.g., `pixelFlags_thing`, not `base_PixelFlags_flag_thing`
        self.applyContext(DrpContext)


class NumStreakDiaSourcesMetric(AnalysisTool):
    """Count DiaSources that fall in a STREAK flag footprint region."""

    def setDefaults(self):
        super().setDefaults()

        # First, select "good" DiaSources that are not obvious garbage
        self.prep.selectors.goodDiaSourceSelector = GoodDiaSourceSelector()
        # Second, select DiaSources with STREAK flag set in the footprint
        self.prep.selectors.flags = FlagSelector(selectWhenTrue=["pixelFlags_streak"])
        # Count them
        self.process.calculateActions.numStreakDiaSources = CountAction(vectorKey="parentDiaSourceId")
        # Set the units for the resulting astropy quantity (ct is count)
        self.produce.metric.units = {"numStreakDiaSources": "ct"}
        # Use, e.g., `pixelFlags_thing`, not `base_PixelFlags_flag_thing`
        self.applyContext(DrpContext)


class NumStreakCenterDiaSourcesMetric(AnalysisTool):
    """Count DiaSources that have the STREAK flag in the center of the source."""

    def setDefaults(self):
        super().setDefaults()

        # First, select "good" DiaSources that are not obvious garbage
        self.prep.selectors.goodDiaSourceSelector = GoodDiaSourceSelector()
        # Second, select DiaSources with STREAK flag set in the source center
        self.prep.selectors.flags = FlagSelector(selectWhenTrue=["pixelFlags_streakCenter"])
        # Count them
        self.process.calculateActions.numStreakCenterDiaSources = CountAction(vectorKey="parentDiaSourceId")
        # Set the units for the resulting astropy quantity (ct is count)
        self.produce.metric.units = {"numStreakCenterDiaSources": "ct"}
        # Use, e.g., `pixelFlags_thing`, not `base_PixelFlags_flag_thing`
        self.applyContext(DrpContext)
