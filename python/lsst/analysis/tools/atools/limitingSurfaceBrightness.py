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

__all__ = ("LimitingSurfaceBrightnessHistPlot",)

from ..actions.plot.histPlot import HistPanel, HistPlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class LimitingSurfaceBrightnessHistPlot(AnalysisTool):
    # We want these to be done and named by band
    parameterizedBand: bool = True

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.hist_limiting_surface_brightness = LoadVector(
            vectorKey="limiting_surface_brightness"
        )

    def finalize(self):
        self.produce.plot = HistPlot()

        self.produce.plot.panels["panel_flux"] = HistPanel()
        self.produce.plot.panels["panel_flux"].label = r"$\mu_{\rm lim}$ (ABmag arcsec$^{-2}$)"
        self.produce.plot.panels["panel_flux"].hists = dict(
            hist_limiting_surface_brightness=r"$\mu_{\rm lim}$",
        )
        self.produce.plot.panels["panel_flux"].rangeType = "sigmaMad"
        self.produce.plot.panels["panel_flux"].lowerRange = 3.5
        self.produce.plot.panels["panel_flux"].upperRange = 3.5
        self.produce.plot.panels["panel_flux"].validate()
        super().finalize()
