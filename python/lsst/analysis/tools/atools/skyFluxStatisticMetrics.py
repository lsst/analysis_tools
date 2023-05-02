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

from lsst.pex.config import Field

from ..actions.scalar.scalarActions import MeanAction, MedianAction, SigmaMadAction, StdevAction
from ..actions.vector.selectors import SkyObjectSelector, SkySourceSelector
from ..interfaces import AnalysisTool
from .coaddVisit import CoaddVisitConfig


class SkyFluxStatisticMetric(AnalysisTool, CoaddVisitConfig):
    """Calculate sky flux statistics.

    This uses the 9-pixel aperture flux for sky sources/objects, and returns
    multiple statistics on the measured fluxes. Note that self.context must be
    set to "visit" (for measurement on sourceTables) or "coadd"
    (for measurement on objectTables).
    """

    fluxType = Field[str](doc="Key to use to retrieve flux column", default="ap09Flux")

    def _setActions(self) -> None:
        if self.context == "coadd":
            name = f"{{band}}_{self.fluxType}"
            self.prep.selectors.skyObjectSelector = SkyObjectSelector()
            self.prep.selectors.skyObjectSelector.bands = []

            # Need to pass a mapping of new names so the default names get the
            # band prepended. Otherwise, each subsequent band's metric will
            # overwrite the current one (e.g., running with g, r bands without
            # this, you would get "meanSky," "meanSky"; with it: "g_meanSky,"
            # "r_meanSky").
            self.produce.metric.newNames = {
                "medianSky": "{band}_medianSky",
                "meanSky": "{band}_meanSky",
                "stdevSky": "{band}_stdevSky",
                "sigmaMADSky": "{band}_sigmaMADSky",
            }
        elif self.context == "visit":
            name = f"{self.fluxType}"
            self.prep.selectors.skySourceSelector = SkySourceSelector()
        else:
            raise ValueError(f"Unsupported {self.context=}")

        self.process.calculateActions.medianSky = MedianAction(vectorKey=name)
        self.process.calculateActions.meanSky = MeanAction(vectorKey=name)
        self.process.calculateActions.stdevSky = StdevAction(vectorKey=name)
        self.process.calculateActions.sigmaMADSky = SigmaMadAction(vectorKey=name)

    def finalize(self) -> None:
        AnalysisTool(self).finalize()
        if not hasattr(self.process.calculateActions, "medianSky"):
            self._setActions()

    def setDefaults(self):
        super().setDefaults()

        self.produce.metric.units = {
            "medianSky": "nJy",
            "meanSky": "nJy",
            "stdevSky": "nJy",
            "sigmaMADSky": "nJy",
        }
