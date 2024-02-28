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
    "CalibStatisticFocalPlanePlot",
    "PtcGainFP",
    "PtcNoiseFP",
    "PtcA00FP",
    "PtcTurnoffFP",
    "PtcMaxRawMeansFP",
    "PtcRowMeanVarianceSlopeFP",
)

from lsst.pex.config import Field

from ..actions.plot.focalPlanePlot import FocalPlaneGeometryPlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class CalibrationTool(AnalysisTool):
    """Class to generate common calibration metrics for value/scatter
    quantities.
    """

    parameterizedBand: bool = False

    def setDefaults(self):
        self.process.buildActions.x = LoadVector(vectorKey="detector")
        self.process.buildActions.y = LoadVector(vectorKey="amplifier")
        self.process.buildActions.detector = LoadVector(vectorKey="detector")
        self.process.buildActions.amplifier = LoadVector(vectorKey="amplifier")
        self.process.buildActions.z = LoadVector()

        self.produce.plot = FocalPlaneGeometryPlot()
        self.produce.plot.statistic = "median"


class CalibStatisticFocalPlanePlot(CalibrationTool):
    """Generates a plot of the focal plane, color-coded according to the
    median of a given measurement (default: "biasMean") on a per-amp basis.
    The median is across multiple bias exposures.
    """

    quantityKey = Field[str](
        default="biasMean", doc="VectorKey to perform the statistic on and to plot per amp and per detector."
    )
    unit = Field[str](default="ADU", doc="Unit of quantity for including on z-axis label.")

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.z.vectorKey = "biasMean"

        self.produce.plot.statistic = "median"
        self.produce.plot.zAxisLabel = "Median of biasMean"

    def finalize(self):
        self.process.buildActions.z.vectorKey = self.quantityKey
        zAxislabel = f"{self.produce.plot.statistic} of {self.quantityKey} ({self.unit})"
        self.produce.plot.zAxisLabel = zAxislabel.capitalize()


class PtcGainFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.z.vectorKey = "ptcGain"
        self.produce.plot.zAxisLabel = "PTC Gain (e-/ADU)"
        self.produce.metric.newNames = {"z": "PTC_GAIN"}


class PtcNoiseFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.z.vectorKey = "ptcNoise"
        self.produce.plot.zAxisLabel = "PTC Readout Noise (ADU^2)"
        self.produce.metric.newNames = {"z": "PTC_NOISE"}


class PtcA00FP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.z.vectorKey = "ptcBfeA00"
        self.produce.plot.zAxisLabel = "PTC BFE A00 (1/e-)"
        self.produce.metric.newNames = {"z": "PTC_BFE_A00"}


class PtcTurnoffFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.z.vectorKey = "ptcTurnoff"
        self.produce.plot.zAxisLabel = "PTC turnoff (ADU)"
        self.produce.metric.newNames = {"z": "PTC_TURNOFF"}


class PtcMaxRawMeansFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.z.vectorKey = "ptcMaxRawMeans"
        self.produce.plot.zAxisLabel = "PTC Maximum of Raw Mean Flux (ADU)"
        self.produce.metric.newNames = {"z": "PTC_MAX_RAW_MEANS"}


class PtcRowMeanVarianceSlopeFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.z.vectorKey = "ptcRowMeanVarianceSlope"
        self.produce.plot.zAxisLabel = "PTC slope of row means vs variance (e-)"
        self.produce.metric.newNames = {"z": "PTC_ROW_MEAN_VARIANCE_SLOPE"}
