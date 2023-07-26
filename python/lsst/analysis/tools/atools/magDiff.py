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

__all__ = ("CModelSubPsfMagMeasSkyGalaxies",)

from ..actions.plot.skyPlot import SkyPlot
from ..actions.vector import CoaddPlotFlagSelector, GalaxySelector, LoadVector, MagDiff, SnSelector
from ..interfaces import AnalysisTool


class CModelSubPsfMagMeasSkyGalaxies(AnalysisTool):
    """Make a plot showing the difference between the cmodel magnitude and the
    PSF magnitude. This plot shows the on sky distribution of these values
    for galaxies.
    """

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["g", "r", "i"]

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.bands = ["i"]
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 100

        self.prep.selectors.galaxySelector = GalaxySelector()
        self.prep.selectors.galaxySelector.vectorKey = "{band}_extendedness"

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.xGalaxies = LoadVector()
        self.process.buildActions.xGalaxies.vectorKey = "coord_ra"
        self.process.buildActions.yGalaxies = LoadVector()
        self.process.buildActions.yGalaxies.vectorKey = "coord_dec"
        self.process.buildActions.patch = LoadVector()
        self.process.buildActions.patch.vectorKey = "patch"
        self.process.buildActions.galaxyStatMask = SnSelector()
        self.process.buildActions.galaxyStatMask.fluxType = "{band}_psfFlux"

        self.process.buildActions.zGalaxies = MagDiff()
        self.process.buildActions.zGalaxies.col1 = "{band}_cModelFlux"
        self.process.buildActions.zGalaxies.col2 = "{band}_psfFlux"

        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["galaxies"]
        self.produce.plot.plotName = "CModel_sub_PSF_meas_galaxies_{band}"
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"
        self.produce.plot.zAxisLabel = "CModel - PSF (mmag) (Meas)"
        self.produce.plot.plotOutlines = False
