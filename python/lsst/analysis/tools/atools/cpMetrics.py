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

from ..actions.plot.focalPlanePlot import FocalPlaneGeometryPlot
from ..actions.scalar.scalarActions import MedianAction, SigmaMadAction
from ..actions.vector import LoadVector, VectorSelector
from ..interfaces import AnalysisTool

__all__ = (
    "BiasScalarMetrics",
    "BiasResidualMeanMedianFP",
    "BiasResidualMeanStdevFP",
    "BiasReadNoiseMedianFP",
    "BiasReadNoiseStdevFP",
    "DarkScalarMetrics",
    "DarkResidualMeanMedianFP",
    "DarkResidualMeanStdevFP",
    "DarkReadNoiseMedianFP",
    "DarkReadNoiseStdevFP",
    "FlatScalarMetrics",
    "FlatResidualMeanMedianFP",
    "FlatResidualMeanStdevFP",
    "FlatNoiseMedianFP",
    "FlatNoiseStdevFP",
    "PtcGainFP",
    "PtcNoiseFP",
    "PtcA00FP",
    "PtcTurnoffFP",
)


class CalibrationMetric(AnalysisTool):
    parameterizedBand: bool = False

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        cls.units = {}
        cls.newNames = {}

    def addPair(self, vectorKey, name, longName):
        setattr(self.process.calculateActions, f"{name}Median", MedianAction(vectorKey=vectorKey))
        setattr(self.process.calculateActions, f"{name}Sigma", SigmaMadAction(vectorKey=vectorKey))

        self.units.update({f"{name}Median": "adu", f"{name}Sigma": "adu"})
        self.newNames.update({f"{name}Median": f"{longName}_median", f"{name}Sigma": f"{longName}_sigmaMad"})

    def addFpPlot(self, vectorKey, statistic, label):
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "detector"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "amplifier"
        self.process.buildActions.detector = LoadVector()
        self.process.buildActions.detector.vectorKey = "detector"
        self.process.buildActions.amplifier = LoadVector()
        self.process.buildActions.amplifier.vectorKey = "amplifier"
        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = vectorKey
        self.process.buildActions.statMask = VectorSelector()
        self.process.buildActions.statMask.vectorKey = "statMask"

        self.produce.plot = FocalPlaneGeometryPlot()
        self.produce.plot.statistic = statistic
        self.produce.plot.xAxisLabel = "x (focal plane)"
        self.produce.plot.yAxisLabel = "y (focal plane)"
        self.produce.plot.zAxisLabel = label


class BiasScalarMetrics(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addPair("MEAN", "mean", "bias_residual_mean")
        self.addPair("NOISE", "noise", "bias_residual_noise")
        self.addPair("CR_NOISE", "crnoise", "bias_residual_crnoise")
        self.addPair("READ_NOISE", "rn", "bias_read_noise")

        self.produce.metric.units = self.units
        self.produce.metric.newNames = self.newNames


class BiasResidualMeanMedianFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("MEAN", "median", "median per-amplifier residual mean (ADU)")


class BiasResidualMeanStdevFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("MEAN", "std", "stdev per-amplifier residual mean (ADU)")


class BiasReadNoiseMedianFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("READ_NOISE", "median", "median read noise (ADU)")


class BiasReadNoiseStdevFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("READ_NOISE", "std", "median read noise (ADU)")


class DarkScalarMetrics(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addPair("MEAN", "mean", "dark_residual_mean")
        self.addPair("NOISE", "noise", "dark_residual_noise")
        self.addPair("CR_NOISE", "crnoise", "dark_residual_crnoise")
        self.addPair("READ_NOISE", "rn", "dark_read_noise")

        self.produce.metric.units = self.units
        self.produce.metric.newNames = self.newNames


class DarkResidualMeanMedianFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("MEAN", "median", "median per-amplifier residual mean (ADU)")


class DarkResidualMeanStdevFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("MEAN", "std", "stdev per-amplifier residual mean (ADU)")


class DarkReadNoiseMedianFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("READ_NOISE", "median", "median read noise (ADU)")


class DarkReadNoiseStdevFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("READ_NOISE", "std", "median read noise (ADU)")


class FlatScalarMetrics(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addPair("MEAN", "mean", "flat_residual_mean")
        self.addPair("NOISE", "noise", "flat_residual_noise")

        self.produce.metric.units = self.units
        self.produce.metric.newNames = self.newNames


class FlatResidualMeanMedianFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("MEAN", "median", "median per-amplifier residual mean (ADU)")


class FlatResidualMeanStdevFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("MEAN", "std", "stdev per-amplifier residual mean (ADU)")


class FlatNoiseMedianFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("NOISE", "median", "median noise (ADU)")


class FlatNoiseStdevFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        self.addFpPlot("NOISE", "std", "median noise (ADU)")


class PtcGainFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("PTC_GAIN", "median", "PTC gain")


class PtcNoiseFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("PTC_NOISE", "median", "PTC read noise")


class PtcA00FP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("PTC_BFE_A00", "median", "PTC A00")


class PtcTurnoffFP(CalibrationMetric):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("PTC_TURNOFF", "median", "PTC turnoff")
