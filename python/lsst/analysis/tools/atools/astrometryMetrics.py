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

from ..actions.scalar import (
    FullRangeAction,
    MaxAction,
    MeanAction,
    MedianAction,
    MinAction,
    SigmaMadAction,
    StdevAction,
)
from ..actions.vector import AngularSeparation, DivideVector
from ..interfaces import AnalysisTool

__all__ = ("AstrometryStatistics",)


class AstrometryStatistics(AnalysisTool):
    """Calculate astrometry metrics from the visit_summary table."""

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.cornersep = AngularSeparation(
            raKey_A="raCorners_0",
            decKey_A="decCorners_0",
            raKey_B="raCorners_2",
            decKey_B="decCorners_2",
            outputUnit="arcminute",
        )
        self.process.buildActions.ratio = DivideVector()
        self.process.buildActions.ratio.actionA = AngularSeparation(
            raKey_A="raCorners_0",
            decKey_A="decCorners_0",
            raKey_B="raCorners_2",
            decKey_B="decCorners_2",
            outputUnit="arcminute",
        )
        self.process.buildActions.ratio.actionB = AngularSeparation(
            raKey_A="raCorners_1",
            decKey_A="decCorners_1",
            raKey_B="raCorners_3",
            decKey_B="decCorners_3",
            outputUnit="arcminute",
        )

        self.process.calculateActions.minCornerSeparation = MinAction(vectorKey="cornersep")
        self.process.calculateActions.maxCornerSeparation = MaxAction(vectorKey="cornersep")
        self.process.calculateActions.minCornerSeparationRatio = MinAction(vectorKey="ratio")
        self.process.calculateActions.maxCornerSeparationRatio = MaxAction(vectorKey="ratio")
        self.process.calculateActions.minPixelScale = MinAction(vectorKey="pixelScale")
        self.process.calculateActions.maxPixelScale = MaxAction(vectorKey="pixelScale")
        self.process.calculateActions.fullRangePixelScale = FullRangeAction(vectorKey="pixelScale")
        self.process.calculateActions.medianPixelScale = MedianAction(vectorKey="pixelScale")
        self.process.calculateActions.sigmaMADPixelScale = SigmaMadAction(vectorKey="pixelScale")
        self.process.calculateActions.meanPixelScale = MeanAction(vectorKey="pixelScale")
        self.process.calculateActions.stdevPixelScale = StdevAction(vectorKey="pixelScale")

        self.produce.metric.units = {
            "minCornerSeparation": "arcmin",
            "maxCornerSeparation": "arcmin",
            "minCornerSeparationRatio": "",
            "maxCornerSeparationRatio": "",
            "minPixelScale": "arcsec",
            "maxPixelScale": "arcsec",
            "medianPixelScale": "arcsec",
            "sigmaMADPixelScale": "arcsec",
            "meanPixelScale": "arcsec",
            "stdevPixelScale": "arcsec",
        }
