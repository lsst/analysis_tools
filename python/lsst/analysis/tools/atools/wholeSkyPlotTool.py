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

__all__ = ("WholeSkyPlotTool",)


from lsst.pex.config import Field

from ..actions.keyedData import KeyedDataUnitAccessAction
from ..actions.plot.wholeSkyPlot import WholeSkyPlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class WholeSkyPlotTool(AnalysisTool):
    """
    WholeSkyPlot of per-tract metric values.

    x and y axes are Right Ascension and Declination, respectively.
    The z-axis colorbar indicates the metric value.
    """

    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter?", default=True
    )

    metric = Field[str](
        doc="The name of the metric to plot.",
        default="e1Diff_{band}_highSNStars_median",
    )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.tract = LoadVector()
        self.process.buildActions.tract.vectorKey = "tract"
        self.produce.plot = WholeSkyPlot()

    def finalize(self):
        self.process.buildActions.z = LoadVector(vectorKey=self.metric)
        self.process.buildActions.zUnit = KeyedDataUnitAccessAction(key=self.metric)
        self.produce.plot.zAxisLabel = self.metric
        sequentialMetrics = ["count", "ean", "edian", "num", "igma", "tdev", "Repeat"]
        if any(sequentialMetric in self.metric for sequentialMetric in sequentialMetrics):
            self.produce.plot.colorMapType = "sequential"
        super().finalize()
