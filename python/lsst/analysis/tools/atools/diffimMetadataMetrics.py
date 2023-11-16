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
    "NUnmergedDiaSourcesMetric",
    "NMergedDiaSourcesMetric",
    "NGoodPixelsMetric",
    "NBadPixelsMetric",
    "NPixelsDetectedPositiveMetric",
    "NPixelsDetectedNegativeMetric",
    "NBadPixelsDetectedPositiveMetric",
    "NBadPixelsDetectedNegativeMetric",
    "SciencePsfSizeMetric",
    "TemplatePsfSizeMetric",
    "ScienceVarianceScaleMetric",
    "TemplateVarianceScaleMetric",
    "TemplateCoverageMetric",
)

from ..actions.scalar import ValueAction
from ..interfaces import AnalysisTool


class NUnmergedDiaSourcesMetric(AnalysisTool):
    """The raw number of DIA source footprints."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nUnmergedDiaSources = ValueAction(vectorKey="nUnmergedDiaSources")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nUnmergedDiaSources": "ct"}


class NMergedDiaSourcesMetric(AnalysisTool):
    """The number of DIA Sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nMergedDiaSources = ValueAction(vectorKey="nMergedDiaSources")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nMergedDiaSources": "ct"}


class NGoodPixelsMetric(AnalysisTool):
    """The number of unflagged pixels in the difference image."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nGoodPixels = ValueAction(vectorKey="nGoodPixels")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nGoodPixels": "ct"}


class NBadPixelsMetric(AnalysisTool):
    """The number of flagged pixels in the difference image."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nBadPixels = ValueAction(vectorKey="nBadPixels")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nBadPixels": "ct"}


class NPixelsDetectedPositiveMetric(AnalysisTool):
    """The number of pixels in the footprints of positive sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nPixelsDetectedPositive = ValueAction(
            vectorKey="nPixelsDetectedPositive"
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nPixelsDetectedPositive": "ct"}


class NPixelsDetectedNegativeMetric(AnalysisTool):
    """The number of pixels in the footprints of negative sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nPixelsDetectedNegative = ValueAction(
            vectorKey="nPixelsDetectedNegative"
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nPixelsDetectedNegative": "ct"}


class NBadPixelsDetectedPositiveMetric(AnalysisTool):
    """The number of flagged pixels in the footprints of positive sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nBadPixelsDetectedPositive = ValueAction(
            vectorKey="nBadPixelsDetectedPositive"
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nBadPixelsDetectedPositive": "ct"}


class NBadPixelsDetectedNegativeMetric(AnalysisTool):
    """The number of flagged pixels in the footprints of negative sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.nBadPixelsDetectedNegative = ValueAction(
            vectorKey="nBadPixelsDetectedNegative"
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"nBadPixelsDetectedNegative": "ct"}


class SciencePsfSizeMetric(AnalysisTool):
    """PSF size of the science image."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.sciencePsfSize = ValueAction(vectorKey="sciencePsfSize")

        # the units for the quantity
        self.produce.metric.units = {"sciencePsfSize": "pixel"}


class TemplatePsfSizeMetric(AnalysisTool):
    """PSF size of the template image."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.templatePsfSize = ValueAction(vectorKey="templatePsfSize")

        # the units for the quantity
        self.produce.metric.units = {"templatePsfSize": "pixel"}


class ScienceVarianceScaleMetric(AnalysisTool):
    """Factor from ScaleVarianceTask for the science image."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.scaleScienceVariance = ValueAction(
            vectorKey="scaleScienceVarianceFactor"
        )

        # the units for the quantity
        self.produce.metric.units = {"scaleScienceVariance": "pixel"}


class TemplateVarianceScaleMetric(AnalysisTool):
    """Factor from ScaleVarianceTask for the template image."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.scaleTemplateVariance = ValueAction(
            vectorKey="scaleTemplateVarianceFactor"
        )

        # the units for the quantity
        self.produce.metric.units = {"scaleTemplateVariance": "pixel"}


class TemplateCoverageMetric(AnalysisTool):
    """Percent of pixels with data in the template image."""

    def setDefaults(self):
        super().setDefaults()

        # Count the number of dia sources
        self.process.calculateActions.templateCoverage = ValueAction(vectorKey="TemplateCoveragePercent")

        # the units for the quantity
        self.produce.metric.units = {"templateCoverage": "percent"}
