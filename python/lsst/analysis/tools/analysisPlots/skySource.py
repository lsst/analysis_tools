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

__all__ = ("SkySourceSkyPlot", "SkySourceHistPlot")

from ..actions.plot.histPlot import HistPanel, HistPlot
from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import LoadVector, SNCalculator
from ..actions.vector.selectors import SkySourceSelector, SnSelector
from ..interfaces import AnalysisPlot


class SkySourceSkyPlot(AnalysisPlot):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.skySourceSelector = SkySourceSelector()

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "coord_ra"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "coord_dec"
        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "ap09Flux"
        self.process.buildActions.statMask = SnSelector()
        self.process.buildActions.statMask.threshold = -1e12
        self.process.buildActions.statMask.fluxType = "psfFlux"

        self.produce = SkyPlot()
        self.produce.plotTypes = ["any"]
        self.produce.plotName = "skySource"
        self.produce.xAxisLabel = "R.A. (degrees)"
        self.produce.yAxisLabel = "Dec. (degrees)"
        self.produce.zAxisLabel = "Sky Source ap09Flux (nJy)"
        self.produce.plotOutlines = False
        self.produce.fixAroundZero = True


class SkySourceHistPlot(AnalysisPlot):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.skySourceSelector = SkySourceSelector()

        self.process.buildActions.hist_psf_flux = LoadVector(vectorKey="psfFlux")
        self.process.buildActions.hist_09_flux = LoadVector(vectorKey="ap09Flux")
        self.process.buildActions.hist_psf_sn = SNCalculator(fluxType="psfFlux")
        self.process.buildActions.hist_09_sn = SNCalculator(fluxType="ap09Flux")

        self.produce = HistPlot()

        self.produce.panels["panel_flux"] = HistPanel()
        self.produce.panels["panel_flux"].label = "Flux (nJy)"
        self.produce.panels["panel_flux"].hists = dict(hist_psf_flux="psfFlux", hist_09_flux="ap09Flux")

        self.produce.panels["panel_sn"] = HistPanel()
        self.produce.panels["panel_sn"].label = "S/N"
        self.produce.panels["panel_sn"].hists = dict(hist_psf_sn="psfFlux S/N", hist_09_sn="ap09Flux S/N")
