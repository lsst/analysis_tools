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

__all__ = (
    "NumFoundFakesDiaAllMetric",
    "FractionFoundFakesDiaAllMetric",
    "NumFakesInjectedAllMetric"
)

from ..actions.scalar import CountAction, FracThreshold
from ..actions.vector import ThresholdSelector
from ..interfaces import AnalysisTool


class NumFoundFakesDiaAllMetric(AnalysisTool):
    """Calculate the number of Fake sources found in the DIA Sources."""

    def setDefaults(self):
        super().setDefaults()

        # filter the fakes that were actually found
        self.prep.selectors.foundFakesSelector = ThresholdSelector(
            vectorKey="diaSourceId",
            op="gt",
            threshold=0)
        # Count the number of fake sources found in the images
        self.process.calculateActions.numFoundFakesDiaSourcesAll = CountAction(vectorKey="diaSourceId")
        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"numFoundFakesDiaSourcesAll": "ct"}


class FractionFoundFakesDiaAllMetric(AnalysisTool):
    """Calculate the number of Fake sources found in the DIA Sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of fake sources found in the images
        self.process.calculateActions.fractionFoundFakesDiaAll = FracThreshold(
            vectorKey="diaSourceId",
            op='gt',
            threshold=0)
        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"fractionFoundFakesDiaAll": "ct"}


class NumFakesInjectedAllMetric(AnalysisTool):
    """Calculate the number of Fake sources found in the DIA Sources."""

    def setDefaults(self):
        super().setDefaults()

        # count the number of total fake sources injected in the images
        self.process.calculateActions.numFakesInjectedAll = CountAction(vectorKey='fakeId')
        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {'numFakesInjectedAll': 'ct'}
