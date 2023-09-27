# This file is part of analysis_ap.
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

__all__ = ("SimpleDiaPlot",)

from ..actions.plot.diaSkyPlot import DiaSkyPanel, DiaSkyPlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class SimpleDiaPlot(AnalysisTool):
    """Single panel DiaSkyPlot for plotting RA/Dec of DiaSources on the sky."""

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.ras = LoadVector()
        self.process.buildActions.ras.vectorKey = "ra"
        self.process.buildActions.decs = LoadVector()
        self.process.buildActions.decs.vectorKey = "dec"
        # TODO: update column name to 'dec' once column names are standardized,
        # i.e., RFC-863

        self.produce.plot = DiaSkyPlot()

        self.produce.plot.panels["panel_main"] = DiaSkyPanel()
        self.produce.plot.panels["panel_main"].xlabel = "RA (deg)"
        self.produce.plot.panels["panel_main"].ylabel = "Dec (deg)"
        self.produce.plot.panels["panel_main"].ra = "ras"
        self.produce.plot.panels["panel_main"].dec = "decs"
        self.produce.plot.panels["panel_main"].rightSpinesVisible = False
