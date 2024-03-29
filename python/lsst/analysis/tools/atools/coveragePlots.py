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

__all__ = ("BaseMultiVisitCoveragePlot", "FocalPlaneMultiVisitCoveragePlot", "RaDecMultiVisitCoveragePlot")

from ..actions.plot.multiVisitCoveragePlot import MultiVisitCoveragePlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class BaseMultiVisitCoveragePlot(AnalysisTool):
    """Base class for multi-visit coverage plots."""

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.produce.plot = MultiVisitCoveragePlot()
        self.updateLoadList()

    # Make sure to call this in a pipeline if the vectorsToLoadList config is
    # modified there.
    def updateLoadList(self):
        for vector in self.produce.plot.vectorsToLoadList:
            setattr(self.process.buildActions, vector, LoadVector())
        for buildAction, vector in zip(self.process.buildActions, self.produce.plot.vectorsToLoadList):
            setattr(buildAction, "vectorKey", vector)


class FocalPlaneMultiVisitCoveragePlot(BaseMultiVisitCoveragePlot):
    """Multi-visit coverage plot in focal plane coordinates."""

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.projection = "focalPlane"
        self.produce.plot.plotName = "focalPlaneMultiVisitCoverage"


class RaDecMultiVisitCoveragePlot(BaseMultiVisitCoveragePlot):
    """Multi-visit coverage plot in RA/Dec coordinates."""

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.plotName = "raDecMultiVisitCoverage"
