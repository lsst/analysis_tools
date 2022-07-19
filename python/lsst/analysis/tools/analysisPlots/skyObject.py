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

__all__ = ("SkyObjectSkyPlot",)

from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import LoadVector
from ..actions.vector.selectors import FlagSelector, SnSelector
from ..interfaces import AnalysisPlot


class SkyObjectSkyPlot(AnalysisPlot):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = FlagSelector()
        self.prep.selectors.flagSelector.selectWhenTrue = ["sky_object"]
        self.prep.selectors.flagSelector.selectWhenFalse = ["{band}_pixelFlags_edge"]

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "coord_ra"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "coord_dec"
        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.threshold = -1e12
        self.process.buildActions.statMask.fluxType = "{band}_psfFlux"

        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "{band}_ap09Flux"

        self.post_process = SkyPlot()
        self.post_process.plotTypes = ["any"]
        self.post_process.plotName = "skyObject_{band}"
        self.post_process.xAxisLabel = "R.A. (degrees)"
        self.post_process.yAxisLabel = "Dec. (degrees)"
        self.post_process.zAxisLabel = "Sky Object Flux (nJy)"
        self.post_process.plotOutlines = False
