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

from lsst.analysis.tools.actions.scalar import MaxAction, MedianAction, MinAction

from ..interfaces import AnalysisMetric

__all__ = ("SeeingStatistics",)


class SeeingStatistics(AnalysisMetric):
    """Calculate seeing statistics from the ccdVisit table.
    Statistics include mean, max, and min.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.scalarActions.medianSeeing = MedianAction(vectorKey=self.vectorKey)
        self.process.scalarActions.minSeeing = MinAction(vectorKey=self.vectorKey)
        self.process.scalarActions.maxSeeing = MaxAction(vectorKey=self.vectorKey)

        self.produce.units = {
            "medianSeeing": "arcsec",
            "minSeeing": "arcsec",
            "maxSeeing": "arcsec",
        }
