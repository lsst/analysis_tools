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

__all__ = ("MatchedRefCoaddTool", "MatchedRefCoaddDiffMagTool")

from ..actions.vector.selectors import GalaxySelector, StarSelector
from ..actions.vector.vectorActions import (
    DivideVector,
    DownselectVector,
    LoadVector,
    MagColumnNanoJansky,
    SubtractVector,
    VectorSelector,
)
from ..interfaces import AnalysisTool


class MatchedRefCoaddTool(AnalysisTool):
    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.fluxes_ref = LoadVector(vectorKey="refcat_flux_{band}")
        # TODO: Why won't vectorKey="fluxes_ref" work?
        # Does it need to be a filterAction?
        self.process.buildActions.mags_ref = MagColumnNanoJansky(
            vectorKey=self.process.buildActions.fluxes_ref.vectorKey
        )

        self.process.buildActions.galaxySelector = GalaxySelector(vectorKey="refExtendedness")
        self.process.buildActions.starSelector = StarSelector(vectorKey="refExtendedness")

        self.process.filterActions.xGalaxies = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="starSelector")
        )


class MatchedRefCoaddDiffMagTool(MatchedRefCoaddTool):
    def matchedRefDiffMagContext(self):
        self.process.buildActions.diff = SubtractVector(
            actionA=MagColumnNanoJansky(vectorKey=self.process.buildActions.fluxes_meas.vectorKey),
            actionB=self.process.buildActions.mags_ref,
        )

    def matchedRefDiffFluxChiContext(self):
        self.process.buildActions.diff = DivideVector(
            actionA=SubtractVector(
                actionA=LoadVector(vectorKey=self.process.buildActions.fluxes_meas.vectorKey),
                actionB=LoadVector(vectorKey=self.process.buildActions.fluxes_ref.vectorKey),
            ),
            actionB=LoadVector(vectorKey=f"{self.process.buildActions.fluxes_meas.vectorKey}Err"),
        )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.fluxes_meas = LoadVector(vectorKey="{band}_cModelFlux")
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="starSelector")
        )
