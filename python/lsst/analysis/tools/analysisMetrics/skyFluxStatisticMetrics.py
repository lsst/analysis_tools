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

__all__ = ("SkyFluxStatisticMetric", "SkyFluxVisitStatisticMetric")

from ..actions.scalar.scalarActions import MeanAction, MedianAction, SigmaMadAction, StdevAction
from ..actions.vector.selectors import SkyObjectSelector, SkySourceSelector
from ..interfaces import AnalysisMetric


class SkyFluxStatisticMetric(AnalysisMetric):
    # Calculate sky flux statistics on objectTable
    parameterizedBand: bool = True
    fluxType: str = "ap09Flux"

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.skyObjectSelector = SkyObjectSelector()
        self.prep.selectors.skyObjectSelector.bands = ["{band}"]

        self.process.calculateActions.medianSky = MedianAction(colKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.meanSky = MeanAction(colKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.stdevSky = StdevAction(colKey=f"{{band}}_{self.fluxType}")
        self.process.calculateActions.sigmaMADSky = SigmaMadAction(colKey=f"{{band}}_{self.fluxType}")

        self.post_process.units = {  # type: ignore
            "medianSky": "nJy",
            "meanSky": "nJy",
            "stdevSky": "nJy",
            "sigmaMADSky": "nJy",
        }


class SkyFluxVisitStatisticMetric(AnalysisMetric):
    # Calculate sky flux statistics on sourceTable (i.e., per visit)
    parameterizedBand: bool = False
    fluxType: str = "ap09Flux"

    def setDefaults(self):
        super().setDefaults()

        self.prep.selectors.skySourceSelector = SkySourceSelector()

        self.process.calculateActions.medianSky = MedianAction(colKey=f"{self.fluxType}")
        self.process.calculateActions.meanSky = MeanAction(colKey=f"{self.fluxType}")
        self.process.calculateActions.stdevSky = StdevAction(colKey=f"{self.fluxType}")
        self.process.calculateActions.sigmaMADSky = SigmaMadAction(colKey=f"{self.fluxType}")

        self.produce.units = {  # type: ignore
            "medianSky": "nJy",
            "meanSky": "nJy",
            "stdevSky": "nJy",
            "sigmaMADSky": "nJy",
        }
