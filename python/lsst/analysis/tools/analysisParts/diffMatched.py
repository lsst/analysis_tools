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

__all__ = ("setMatchedRefCoaddDefaults", "setMatchedRefCoaddDiffMagDefaults")

from ..interfaces import AnalysisTool
from ..actions.vector.selectors import GalaxySelector, StarSelector
from ..actions.vector.vectorActions import (
    DownselectVector,
    MagColumnNanoJansky,
    MagDiff,
    VectorSelector,
)


def setMatchedRefCoaddDefaults(tool: AnalysisTool):
    # tool.prep.selectors.flagSelector = CoaddPlotFlagSelector()
    tool.process.buildActions.mags_ref = MagColumnNanoJansky(vectorKey="refcat_flux_{band}")

    # pre-compute a stellar selector mask so it can be used in the filter
    # actions while only being computed once, alternatively the stellar
    # selector could be calculated and applied twice in the filter stage
    tool.process.buildActions.galaxySelector = GalaxySelector(columnKey='refExtendedness')
    tool.process.buildActions.starSelector = StarSelector(columnKey='refExtendedness')

    tool.process.filterActions.xGalaxies = DownselectVector(
        vectorKey="mags_ref", selector=VectorSelector(vectorKey="galaxySelector")
    )
    tool.process.filterActions.xStars = DownselectVector(
        vectorKey="mags_ref", selector=VectorSelector(vectorKey="starSelector")
    )


def setMatchedRefCoaddDiffMagDefaults(tool: AnalysisTool, column_flux: str = None):
    if column_flux is None:
        column_flux = "{band}_cModelFlux"
    tool.process.buildActions.magDiff = MagDiff(
        col1=column_flux,
        col2=tool.process.buildActions.mags_ref.vectorKey,
    )

    tool.process.filterActions.yGalaxies = DownselectVector(
        vectorKey="magDiff", selector=VectorSelector(vectorKey="galaxySelector")
    )
    tool.process.filterActions.yStars = DownselectVector(
        vectorKey="magDiff", selector=VectorSelector(vectorKey="starSelector")
    )
