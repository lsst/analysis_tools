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

__all__ = ("LimitingSurfBriXYPlot",)

from ..actions.plot.xyPlot import XYPlot
from ..actions.scalar.scalarActions import IqrHistAction, MedianHistAction
from ..actions.vector import ConstantValue, LoadVector
from ..interfaces import AnalysisTool


class LimitingSurfBriXYPlot(AnalysisTool):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.x.vectorKey = "bin_mid"
        self.process.buildActions.y.vectorKey = "hist"
        self.process.buildActions.xerr = ConstantValue(value=0)
        self.process.buildActions.yerr = ConstantValue(value=0)

        self.process.calculateActions.median = MedianHistAction()
        self.process.calculateActions.median.histKey = "hist"
        self.process.calculateActions.median.midKey = "bin_mid"
        self.process.calculateActions.iqr = IqrHistAction()
        self.process.calculateActions.iqr.histKey = "hist"
        self.process.calculateActions.iqr.midKey = "bin_mid"

        self.produce.plot = XYPlot()
        self.produce.plot.xAxisLabel = r"$mu_{\rm lim}$ (ABmag)"
        self.produce.plot.yAxisLabel = "Frequency"
        self.produce.plot.yScale = "linear"
        self.produce.plot.xLine = 0
        self.produce.plot.strKwargs = {
            "fmt": "-",
            "color": "b",
        }

        self.produce.metric.units = {
            "median": "mag",
            "iqr": "mag",
        }

        self.produce.metric.newNames = {
            "median": "limiting_surfBri_median",
            "iqr": "limiting_surfBri_iqr",
        }