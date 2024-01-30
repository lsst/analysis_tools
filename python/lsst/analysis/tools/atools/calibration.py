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

__all__ = ("MedReadNoiseFocalPlanePlot",)

from ..actions.plot.focalPlanePlot import FocalPlaneGeometryPlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class MedReadNoiseFocalPlanePlot(AnalysisTool):
    """Generates a plot of the focal plane, color-coded according to the
    median bias read noise on a per-amp basis. The median is across
    multiple bias exposures.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "biasReadNoise"
        self.process.buildActions.detector = LoadVector()
        self.process.buildActions.detector.vectorKey = "detector"
        self.process.buildActions.amplifier = LoadVector()
        self.process.buildActions.amplifier.vectorKey = "amplifier"

        self.produce.plot = FocalPlaneGeometryPlot()
        # TO DO: In anticipation of the addHistogram option
        # self.produce.plot.addHistogram = True
        self.produce.plot.xAxisLabel = "x (mm)"
        self.produce.plot.yAxisLabel = "y (mm)"
        self.produce.plot.zAxisLabel = "Med. Readnoise"
        self.produce.plot.statistic = "median"
