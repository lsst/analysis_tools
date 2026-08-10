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

__all__ = ("CalibStatisticFocalPlanePlot", "CalibAmpMetric")

from lsst.pex.config import Field, ListField

from ..actions.plot.focalPlanePlot import FocalPlaneGeometryPlot
from ..actions.scalar.scalarActions import MedianAction
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool, KeyedDataAction, NoPlot, Vector


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


class PerAmpKeyedScalars(KeyedDataAction):
    """Repack per-row values into one mapping keyed by amp."""

    outputKey = Field[str](doc="Name the per-amp mapping is stored under.")
    groupKeys = ListField[str](doc="Row keys joined to form each amp label.",
                               default=["detector", "amplifier"])
    valueKey = Field[str](doc="Per-row value to turn into per-amp metrics.")

    def getInputSchema(self):
        for k in self.groupKeys:
            yield (k, Vector)
        yield (self.valueKey, Vector)

    def __call__(self, data, **kwargs):
        value = data[self.valueKey.format(**kwargs)]
        cols = [data[k.format(**kwargs)] for k in self.groupKeys]
        mapping = {"_".join(str(c[i]) for c in cols): float(value[i])
                   for i in range(len(value))}
        return {self.outputKey: mapping}   # wrapped so BaseProcess keeps it grouped


class CalibAmpMetric(CalibrationTool):
    """Doc to add
    """
    quantityKey = Field[str](
        default="DEFECTS_DEFECT_PIXELS", doc="VectorKey to output the per amp value from.")
    unit = Field[str](default="", doc="Unit of value per amp.")

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot = NoPlot()

    def finalize(self):
        self.process.calculateActions.perAmp = PerAmpKeyedScalars(
            outputKey="nDefects", valueKey=self.quantityKey)
        self.produce.metric.units = {"nDefects": self.unit}


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

        self.process.calculateActions.median = MedianAction()
        self.process.calculateActions.median.vectorKey = "biasMean"

        self.produce.plot.statistic = "median"
        self.produce.plot.zAxisLabel = "Median of biasMean"

    def finalize(self):
        self.process.buildActions.z.vectorKey = self.quantityKey
        self.process.calculateActions.median.vectorKey = self.quantityKey
        self.produce.metric.units = {"median": "adu"}
        zAxislabel = f"{self.produce.plot.statistic} of {self.quantityKey} ({self.unit})"
        self.produce.plot.zAxisLabel = zAxislabel.capitalize()
