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

__all__ = ("SkyFluxStatisticMetric",)

from ..actions.scalar.scalarActions import MeanAction, MedianAction, SigmaMadAction, StdevAction
from ..actions.vector.selectors import SkyObjectSelector, SkySourceSelector
from ..interfaces import AnalysisMetric


class SkyFluxStatisticMetric(AnalysisMetric):
    """Calculate sky flux statistics. This uses the 9-pixel aperture flux for
    sky sources/objects, and returns multiple statistics on the measured
    fluxes. Note that either visitContext (measurement on sourceTable) or
    coaddContext (measurement on objectTable) must be specified.
    """

    fluxType: str = "ap09Flux"

    def visitContext(self) -> None:
        self.prep.selectors.skySourceSelector = SkySourceSelector()
        self._setActions(f"{self.fluxType}")

    def coaddContext(self) -> None:
        self.prep.selectors.skyObjectSelector = SkyObjectSelector()
        self.prep.selectors.skyObjectSelector.bands = []
        self._setActions(f"{{band}}_{self.fluxType}")

        # Need to pass a mapping of new names so the default names get the
        # band prepended. Otherwise, each subsequent band's metric will
        # overwrite the current one (e.g., running with g, r bands without
        # this, you would get "meanSky," "meanSky"; with it: "g_meanSky,"
        # "r_meanSky").
        self.produce.newNames = {
            "medianSky": "{band}_medianSky",
            "meanSky": "{band}_meanSky",
            "stdevSky": "{band}_stdevSky",
            "sigmaMADSky": "{band}_sigmaMADSky",
        }

    def _setActions(self, name: str) -> None:
        self.process.calculateActions.medianSky = MedianAction(vectorKey=name)
        self.process.calculateActions.meanSky = MeanAction(vectorKey=name)
        self.process.calculateActions.stdevSky = StdevAction(vectorKey=name)
        self.process.calculateActions.sigmaMADSky = SigmaMadAction(vectorKey=name)

    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "medianSky": "nJy",
            "meanSky": "nJy",
            "stdevSky": "nJy",
            "sigmaMADSky": "nJy",
        }
