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

__all__ = ("DiaObjectPlot",)

from ..actions.plot.calculateRange import MinMax
from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import DownselectVector, LoadVector
from ..actions.vector.selectors import FiniteSelector, ThresholdSelector
from ..interfaces import AnalysisTool


class DiaObjectPlot(AnalysisTool):
    """Make a plot of DiaObjects on the sky."""

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.x = LoadVector(vectorKey="ra")
        self.process.buildActions.y = LoadVector(vectorKey="dec")
        self.process.buildActions.z = LoadVector(vectorKey="nDiaSources")

        # statMask is required for SkyPlot
        # it computes the nanMedian and nanSigmaMad for the z array,
        # with some selector applied
        self.process.buildActions.statMask = FiniteSelector(vectorKey="nDiaSources")

        # only plot diaObjects composed of 5 or fewer diaSources
        self.process.filterActions.z = DownselectVector(vectorKey="z")
        self.process.filterActions.z.selector = ThresholdSelector(
            vectorKey="nDiaSources", op="le", threshold=5
        )
        self.process.filterActions.x = DownselectVector(
            vectorKey="x", selector=self.process.filterActions.z.selector
        )
        self.process.filterActions.y = DownselectVector(
            vectorKey="y", selector=self.process.filterActions.z.selector
        )

        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["any"]
        self.produce.plot.plotName = "DiaObjects with 5 or fewer DiaSources"

        self.produce.plot.xAxisLabel = "R.A. (deg)"
        self.produce.plot.yAxisLabel = "Dec. (deg)"
        self.produce.plot.zAxisLabel = "Number of associated DiaSources"

        self.produce.plot.colorbarRange = MinMax
        self.produce.plot.plotOutlines = True
        self.produce.plot.doBinning = False
        self.produce.plot.alpha = 0.2
        self.produce.plot.scatPtSize = 3
