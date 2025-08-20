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

__all__ = ("WholeTractImageTool", "WholeTractClippedMaskTool", "WholeTractPostageStampTool")

from ..actions.plot.calculateRange import Linear, MinMax
from ..actions.plot.wholeTractImage import WholeTractImage
from ..interfaces import AnalysisTool


class WholeTractImageTool(AnalysisTool):
    """Makes a figure displaying all coadd images covering a whole tract."""

    propagateData = True
    parameterizedBand = False

    def setDefaults(self):
        super().setDefaults()

        self.prep.keysToLoad = ["image"]
        self.produce.plot = WholeTractImage()


class WholeTractClippedMaskTool(AnalysisTool):
    """Makes a figure displaying the pixels flagged as CLIPPED for a whole
    tract."""

    propagateData = True
    parameterizedBand = False

    def setDefaults(self):
        super().setDefaults()

        self.prep.keysToLoad = ["mask"]
        self.produce.plot = WholeTractImage()
        self.produce.plot.component = "mask"
        self.produce.plot.bitmaskPlanes = ["CLIPPED"]
        self.produce.plot.interval = MinMax()
        self.produce.plot.stretch = Linear()


class WholeTractPostageStampTool(AnalysisTool):
    """Makes a postage-stamp figure displaying all coadd images covering a
    whole tract."""

    propagateData = True
    parameterizedBand = False

    def setDefaults(self):
        super().setDefaults()

        self.prep.keysToLoad = ["image"]
        self.produce.plot = WholeTractImage()
        self.produce.plot.displayAsPostageStamp = True
