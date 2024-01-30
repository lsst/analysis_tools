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
    "MedReadNoiseFocalPlanePlot",
    "PtcGainFP",
    "PtcNoiseFP",
    "PtcA00FP",
    "PtcTurnoffFP")

from ..actions.plot.focalPlanePlot import FocalPlaneGeometryPlot
from ..actions.scalar.scalarActions import MedianAction, SigmaMadAction
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class CalibrationTool(AnalysisTool):
    """Class to generate common calibration metrics for value/scatter
    quantities (from python/lsst/analysis/tools/atools/cpMetrics.py
    in DM-40473)
    """

    parameterizedBand: bool = False

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        cls.units = {}
        cls.newNames = {}

    def addPair(self, vectorKey, name, longName):
        """Add a pair of value/scatter metrics from a given catalog key.
        Parameters
        ----------
        vectorKey : `str`
            Name of the catalog key to load.
        name : `str`
            Short name for this metric.
        longName : `str`
            Detailed metric name, including calibration stage and product.
        """
        setattr(self.process.calculateActions, f"{name}Median", MedianAction(vectorKey=vectorKey))
        setattr(self.process.calculateActions, f"{name}Sigma", SigmaMadAction(vectorKey=vectorKey))

        self.units.update({f"{name}Median": "adu", f"{name}Sigma": "adu"})
        self.newNames.update({f"{name}Median": f"{longName}_median", f"{name}Sigma": f"{longName}_sigmaMad"})

    def addFpPlot(self, vectorKey, statistic, label):
        """Add focal plan geometry plot.
        Parameters
        ----------
        vectorKey : `str`
             Name of the catalog key to load.
        statistic : `str`
             Statistic to use in binning per-amplifier data points.
        label : `str`
             Label to apply to the output z-axis.
        """
        self.process.buildActions.x = LoadVector(vectorKey="detector")
        self.process.buildActions.y = LoadVector(vectorKey="amplifier")
        self.process.buildActions.detector = LoadVector(vectorKey="detector")
        self.process.buildActions.amplifier = LoadVector(vectorKey="amplifier")
        self.process.buildActions.z = LoadVector(vectorKey=vectorKey)

        self.produce.plot = FocalPlaneGeometryPlot()
        self.produce.plot.statistic = statistic
        self.produce.plot.xAxisLabel = "x (focal plane)"
        self.produce.plot.yAxisLabel = "y (focal plane)"
        self.produce.plot.zAxisLabel = label


class MedReadNoiseFocalPlanePlot(AnalysisTool):
    """Generates a plot of the focal plane, color-coded according to the
    median bias read noise on a per-amp basis. The median is across
    multiple bias exposures.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.z = LoadVector()
        self.process.buildActions.z.vectorKey = "biasReadNoise"
        self.process.buildActions.detector = LoadVector()
        self.process.buildActions.detector.vectorKey = "detector"
        self.process.buildActions.amplifier = LoadVector()
        self.process.buildActions.amplifier.vectorKey = "amplifier"

        self.produce.plot = FocalPlaneGeometryPlot()
        # TO DO: In anticipation of the addHistogram option
        # self.produce.plot.addHistogram = True
        self.produce.plot.xAxisLabel = "x (mm)"
        self.produce.plot.yAxisLabel = "y (mm)"
        self.produce.plot.zAxisLabel = "Med. Readnoise"
        self.produce.plot.statistic = "median"


class PtcGainFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("ptcGain", "median", "PTC gain")


class PtcNoiseFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("ptcNoise", "median", "PTC read noise")


class PtcA00FP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("ptcBfeA00", "median", "PTC BFE A00")


class PtcTurnoffFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("ptcTurnoff", "median", "PTC turnoff")


class PtcMaxRawMeansFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("ptcMaxRawMeans", "median", "PTC max of raw means")


class PtcRowMeanVarianceSlopeFP(CalibrationTool):
    def setDefaults(self):
        super().setDefaults()
        # This should only have one entry, so the statistic doesn't
        # matter much.
        self.addFpPlot("ptcRowMeanVarianceSlope", "median", "PTC slope of row means vs variance")
