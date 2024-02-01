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

__all__ = ("CalibStatisticFocalPlanePlot", "CalibStatisticFocalPlaneMetric")

from lsst.pex.config import Field

from ..actions.vector import LoadVector, MultiCriteriaDownselectVector, ThresholdSelector, BandSelector
from ..actions.scalar import MedianAction
from ..actions.plot import FocalPlaneGeometryPlot
from ..interfaces import AnalysisTool


class CalibStatisticFocalPlanePlot(AnalysisTool):
    """Generates a plot of the focal plane, color-coded according to the
    median of a given measurement (default: "biasMean") on a per-amp basis.
    The median is across multiple bias exposures.
    """

    quantity = Field[str](
        doc="Qauntity on which to perform statistics.",
        default="biasNoise",
    )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "biasMean"
        self.process.buildActions.detector = LoadVector()
        self.process.buildActions.detector.vectorKey = "detector"
        self.process.buildActions.amplifier = LoadVector()
        self.process.buildActions.amplifier.vectorKey = "amplifier"

        self.produce.plot = FocalPlaneGeometryPlot()
        # TO DO: In anticipation of the addHistogram option
        # self.produce.plot.addHistogram = True
        self.produce.plot.statistic = "median"
        self.produce.plot.xAxisLabel = "x (mm)"
        self.produce.plot.yAxisLabel = "y (mm)"
        self.produce.plot.zAxisLabel = "Median of Mean Bias"
        self.produce.plot.statistic = "median"

    def finalize(self):
        self.process.buildActions.z.vectorKey = self.quantity
        zAxislabel = f"{self.produce.plot.statistic} of {self.process.buildActions.z.vectorKey}"
        self.produce.plot.zAxisLabel = zAxislabel.capitalize()


class CalibStatisticFocalPlaneMetric(AnalysisTool):
    """Calculates the median value of a quantity (default: biasNoise) across
    multiple exposures on a per-amp-per-detector basis.
    """

    quantity = Field[str](
        doc="Qauntity on which to perform statistics.",
        default="biasNoise",
    )

    _n_detectors = 1
    _n_amplifiers = 16

    def setDefaults(self):
        super().setDefaults()

        for detectorId in range(self._n_detectors):
            for ampId in range(self._n_amplifiers):
                ampName = f"C{ampId:02}"
                self.process.filterActions.quantity = MultiCriteriaDownselectVector(
                    vectorKey="biasMean",
                )
                self.process.filterActions.quantity.selectors.amp = BandSelector(
                    vectorKey="amplifier",
                    bands=[ampName],
                )
                self.process.filterActions.quantity.selectors.detector = ThresholdSelector(
                    op="eq",
                    threshold=detectorId,
                    vectorKey="detector",
                )

                attrName = f"D{detectorId:03}_C{ampId:02}"
                setattr(
                    self.process.calculateActions,
                    attrName,
                    MedianAction(vectorKey="quantity"),
                )
                self.produce.metric.units[attrName] = "count"

    def finalize(self):
        self.process.filterActions.quantity.vectorKey = self.quantity
