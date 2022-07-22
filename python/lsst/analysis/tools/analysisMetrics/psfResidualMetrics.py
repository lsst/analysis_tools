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

__all__ = ("ShapeSizeFractionalMetric",)

from ..actions.vector import DownselectVector, VectorSelector
from ..analysisParts.shapeSizeFractional import BasePsfResidualMixin
from ..interfaces import AnalysisMetric


class BasePsfResidualMetric(AnalysisMetric, BasePsfResidualMixin):
    def setDefaults(self):
        super().setDefaults()
        self.produce.units = {
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        }


class ShapeSizeFractionalMetric(BasePsfResidualMetric):
    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="fracDiff", selector=VectorSelector(vectorKey="starSelector")
        )
