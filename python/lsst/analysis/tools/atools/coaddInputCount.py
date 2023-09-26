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

__all__ = ("CoaddInputCount", )

from ..actions.plot.skyPlot import SkyPlot
from ..actions.scalar.scalarActions import MeanAction, MedianAction, SigmaMadAction
from ..actions.vector import (
    CoaddPlotFlagSelector,
    LoadVector,
    SnSelector,
)
from ..interfaces import AnalysisTool


class CoaddInputCount(AnalysisTool):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        # Set this to an empty list to look at the band
        # the plot is being made in.
        self.prep.selectors.flagSelector.bands = []

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 100

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "coord_ra"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "coord_dec"
        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.fluxType = "{band}_psfFlux"
        self.process.buildActions.statMask.threshold = 100

        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "{band}_inputCount"

        #self.process.buildActions.patch = LoadVector()
        #self.process.buildActions.patch.vectorKey = "patch"

        self.process.calculateActions.median = MedianAction()
        self.process.calculateActions.median.vectorKey = "z"

        self.process.calculateActions.mean = MeanAction()
        self.process.calculateActions.mean.vectorKey = "z"

        self.process.calculateActions.sigmaMad = SigmaMadAction()
        self.process.calculateActions.sigmaMad.vectorKey = "z"

        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["any"]
        self.produce.plot.plotName = "{band}_inputCount"
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"
        self.produce.plot.zAxisLabel = "Input Count"
        self.produce.plot.plotOutlines = True
        self.produce.plot.addExtremeScatter = False

        self.produce.metric.units = {"median": "ct", "sigmaMad": "ct", "mean": "ct"}

        self.produce.metric.newNames = {
            "median": "{band}_inputCount_median",
            "mean": "{band}_inputCount_mean",
            "sigmaMad": "{band}_inputCount_sigmaMad",
        }
