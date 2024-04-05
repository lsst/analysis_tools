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
__all__ = ("DiffimSpatialMetricsHistPlot", "DiffimSpatialMetricsInterpolatePlot")

from lsst.pex.config import Field

from ..actions.plot.histPlot import HistPanel, HistPlot
from ..actions.plot.interpolateDetectorPlot import InterpolateDetectorMetricPlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class DiffimSpatialMetricsHistPlot(AnalysisTool):
    """Create histograms of the fraction of pixels with certain mask planes
    set.
    """

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.bad_mask_fraction = LoadVector(vectorKey="bad_mask_fraction")
        self.process.buildActions.cr_mask_fraction = LoadVector(vectorKey="cr_mask_fraction")
        self.process.buildActions.detected_mask_fraction = LoadVector(vectorKey="detected_mask_fraction")
        self.process.buildActions.detected_negative_mask_fraction = LoadVector(
            vectorKey="detected_negative_mask_fraction"
        )
        self.process.buildActions.intrp_mask_fraction = LoadVector(vectorKey="intrp_mask_fraction")
        self.process.buildActions.no_data_mask_fraction = LoadVector(vectorKey="no_data_mask_fraction")
        self.process.buildActions.sat_mask_fraction = LoadVector(vectorKey="sat_mask_fraction")
        self.process.buildActions.sat_template_mask_fraction = LoadVector(
            vectorKey="sat_template_mask_fraction"
        )
        self.process.buildActions.streak_mask_fraction = LoadVector(vectorKey="streak_mask_fraction")

        self.produce.plot = HistPlot()

        self.produce.plot.panels["panel_flags"] = HistPanel()
        self.produce.plot.panels["panel_flags"].label = "Flagged pixel fraction"
        self.produce.plot.panels["panel_flags"].bins = 20
        self.produce.plot.panels["panel_flags"].rangeType = "fixed"
        self.produce.plot.panels["panel_flags"].lowerRange = 0
        self.produce.plot.panels["panel_flags"].upperRange = 1.0
        self.produce.plot.panels["panel_flags"].hists = dict(
            no_data_mask_fraction="No data",
            sat_mask_fraction="Saturated",
            bad_mask_fraction="Bad",
            cr_mask_fraction="Cosmic ray",
            detected_mask_fraction="Detected",
            detected_negative_mask_fraction="Detected negative",
            intrp_mask_fraction="Interpolated",
            sat_template_mask_fraction="Saturated template",
            streak_mask_fraction="Streak",
        )


class DiffimSpatialMetricsInterpolatePlot(AnalysisTool):
    """Interpolate spatially-sampled metric values and create low-resolution
    images of the result.
    """

    metricName = Field[str](doc="Metric name to interpolate", optional=False)
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.produce.plot = InterpolateDetectorMetricPlot()
        self.process.buildActions.x = LoadVector()
        self.process.buildActions.x.vectorKey = "x"
        self.process.buildActions.y = LoadVector()
        self.process.buildActions.y.vectorKey = "y"

    def finalize(self):
        self.process.buildActions.metricValues = LoadVector()
        self.process.buildActions.metricValues.vectorKey = self.metricName
