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

__all__ = ("MatchedRefCoaddTool", "MatchedRefCoaddDiffMagTool", "MatchedRefCoaddDiffPositionTool")

from lsst.pex.config import ChoiceField

from ..actions.vector.selectors import GalaxySelector, StarSelector
from ..actions.vector.vectorActions import (
    ConstantValue,
    DivideVector,
    DownselectVector,
    LoadVector,
    MagColumnNanoJansky,
    SubtractVector,
    VectorSelector,
)
from ..interfaces import AnalysisTool


class MatchedRefCoaddTool(AnalysisTool):
    """Base tool for matched-to-reference metrics/plots on coadds.

    Metrics/plots are expected to use the reference magnitude and
    require separate star/galaxy/all source selections.

    Notes
    -----
    The tool does not use a standard coadd flag selector, because
    it is expected that the matcher has been configured to select
    appropriate candidates (and stores a match_candidate column).
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.fluxes_ref = LoadVector(vectorKey="refcat_flux_{band}")
        # TODO: Why won't vectorKey="fluxes_ref" work?
        # Does it need to be a filterAction?
        self.process.buildActions.mags_ref = MagColumnNanoJansky(
            vectorKey=self.process.buildActions.fluxes_ref.vectorKey
        )

        # Select any finite extendedness (but still exclude NaNs)
        self.process.buildActions.allSelector = StarSelector(
            vectorKey="refExtendedness", extendedness_maximum=1.0
        )
        self.process.buildActions.galaxySelector = GalaxySelector(vectorKey="refExtendedness")
        self.process.buildActions.starSelector = StarSelector(vectorKey="refExtendedness")

        self.process.filterActions.xAll = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.xGalaxies = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="starSelector")
        )


class MatchedRefCoaddDiffMagTool(MatchedRefCoaddTool):
    """Base tool for diffs between reference and measured coadd mags.

    Notes
    -----
    The default model flux is cModel.
    """

    def matchedRefDiffContext(self):
        self.process.buildActions.diff = SubtractVector(
            actionA=MagColumnNanoJansky(
                vectorKey=self.process.buildActions.fluxes_meas.vectorKey, returnMillimags=True
            ),
            actionB=DivideVector(
                actionA=self.process.buildActions.mags_ref,
                # To convert to mmag
                actionB=ConstantValue(value=1e-3),
            ),
        )

    def matchedRefChiContext(self):
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
        self.process.filterActions.yAll = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="starSelector")
        )


class MatchedRefCoaddDiffPositionTool(MatchedRefCoaddTool):
    """Base tool for diffs between reference and measured coadd astrometry."""

    variable = ChoiceField[str](
        doc="The astrometric variable to compute metrics for",
        allowed={
            "x": "x",
            "y": "y",
        },
        optional=False,
    )

    # TODO: Determine if this can be put back into setDefaults w/o this:
    # lsst.pex.config.config.FieldValidationError:
    # Field 'process.buildActions.pos_meas.vectorKey' failed validation:
    # Required value cannot be None
    def _setPos(self):
        self.process.buildActions.pos_meas = LoadVector(vectorKey=self.variable)
        self.process.buildActions.pos_ref = LoadVector(vectorKey=f"refcat_{self.variable}")

    def matchedRefDiffContext(self):
        self._setPos()
        self.process.buildActions.diff = SubtractVector(
            actionA=self.process.buildActions.pos_meas,
            actionB=self.process.buildActions.pos_ref,
        )

    def matchedRefChiContext(self):
        self._setPos()
        self.process.buildActions.diff = DivideVector(
            actionA=SubtractVector(
                actionA=self.process.buildActions.pos_meas,
                actionB=self.process.buildActions.pos_ref,
            ),
            actionB=LoadVector(vectorKey=f"{self.process.buildActions.pos_meas.vectorKey}Err"),
        )

    def setDefaults(self):
        super().setDefaults()

        self.process.filterActions.yAll = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="starSelector")
        )
