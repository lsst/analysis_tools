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
    "BiasPercentilePlot",
    "DarkPercentilePlot",
    "FlatPercentilePlot",
)

from ..actions.plot.percentilePlot import PercentilePlot

# from ..actions.scalar.scalarActions import MedianAction, SigmaMadAction
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool
from lsst.pex.config import Field


class BiasPercentilePlot(AnalysisTool):
    """Plot the percentiles of the normalized amplifier bias distributions."""

    parameterizedBand = Field[bool](
        default=False,
        doc="Does this AnalysisTool support band as a name parameter",
    )

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.amplifier = LoadVector()
        self.process.buildActions.amplifier.vectorKey = "amplifier"

        self.process.buildActions.detector = LoadVector()
        self.process.buildActions.detector.vectorKey = "detector"

        self.process.buildActions.percentile_0 = LoadVector()
        self.process.buildActions.percentile_0.vectorKey = "biasDistribution_0.0"

        self.process.buildActions.percentile_5 = LoadVector()
        self.process.buildActions.percentile_5.vectorKey = "biasDistribution_5.0"

        self.process.buildActions.percentile_16 = LoadVector()
        self.process.buildActions.percentile_16.vectorKey = "biasDistribution_16.0"

        self.process.buildActions.percentile_50 = LoadVector()
        self.process.buildActions.percentile_50.vectorKey = "biasDistribution_50.0"

        self.process.buildActions.percentile_84 = LoadVector()
        self.process.buildActions.percentile_84.vectorKey = "biasDistribution_84.0"

        self.process.buildActions.percentile_95 = LoadVector()
        self.process.buildActions.percentile_95.vectorKey = "biasDistribution_95.0"

        self.process.buildActions.percentile_100 = LoadVector()
        self.process.buildActions.percentile_100.vectorKey = "biasDistribution_100.0"

        # self.process.calculateActions.mag50 = Mag50Action()
        # self.process.calculateActions.mag50.vectorKey = "{band}_mag_ref"
        # self.process.calculateActions.mag50.matchDistanceKey = "matchDistance"

        self.produce.plot = PercentilePlot()
        # self.produce.metric.units = {"mag50": "mag"}
        # self.produce.metric.newNames = {"mag50": "{band}_mag50"}


class DarkPercentilePlot(AnalysisTool):
    """Plot the percentiles of the normalized amplifier dark distributions."""

    parameterizedBand = Field[bool](
        default=False,
        doc="Does this AnalysisTool support band as a name parameter",
    )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.amplifier = LoadVector()
        self.process.buildActions.amplifier.vectorKey = "amplifier"

        self.process.buildActions.detector = LoadVector()
        self.process.buildActions.detector.vectorKey = "detector"

        self.process.buildActions.percentile_0 = LoadVector()
        self.process.buildActions.percentile_0.vectorKey = "darkDistribution_0.0"

        self.process.buildActions.percentile_5 = LoadVector()
        self.process.buildActions.percentile_5.vectorKey = "darkDistribution_5.0"

        self.process.buildActions.percentile_16 = LoadVector()
        self.process.buildActions.percentile_16.vectorKey = "darkDistribution_16.0"

        self.process.buildActions.percentile_50 = LoadVector()
        self.process.buildActions.percentile_50.vectorKey = "darkDistribution_50.0"

        self.process.buildActions.percentile_84 = LoadVector()
        self.process.buildActions.percentile_84.vectorKey = "darkDistribution_84.0"

        self.process.buildActions.percentile_95 = LoadVector()
        self.process.buildActions.percentile_95.vectorKey = "darkDistribution_95.0"

        self.process.buildActions.percentile_100 = LoadVector()
        self.process.buildActions.percentile_100.vectorKey = "darkDistribution_100.0"

        # self.process.calculateActions.mag50 = Mag50Action()
        # self.process.calculateActions.mag50.vectorKey = "{band}_mag_ref"
        # self.process.calculateActions.mag50.matchDistanceKey = "matchDistance"

        self.produce.plot = PercentilePlot()
        # self.produce.metric.units = {"mag50": "mag"}
        # self.produce.metric.newNames = {"mag50": "{band}_mag50"}


class FlatPercentilePlot(AnalysisTool):
    """Plot the percentiles of the normalized amplifier flat distributions."""

    parameterizedBand = Field[bool](
        default=False,
        doc="Does this AnalysisTool support band as a name parameter",
    )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.amplifier = LoadVector()
        self.process.buildActions.amplifier.vectorKey = "amplifier"

        self.process.buildActions.detector = LoadVector()
        self.process.buildActions.detector.vectorKey = "detector"

        self.process.buildActions.percentile_0 = LoadVector()
        self.process.buildActions.percentile_0.vectorKey = "flatDistribution_0.0"

        self.process.buildActions.percentile_5 = LoadVector()
        self.process.buildActions.percentile_5.vectorKey = "flatDistribution_5.0"

        self.process.buildActions.percentile_16 = LoadVector()
        self.process.buildActions.percentile_16.vectorKey = "flatDistribution_16.0"

        self.process.buildActions.percentile_50 = LoadVector()
        self.process.buildActions.percentile_50.vectorKey = "flatDistribution_50.0"

        self.process.buildActions.percentile_84 = LoadVector()
        self.process.buildActions.percentile_84.vectorKey = "flatDistribution_84.0"

        self.process.buildActions.percentile_95 = LoadVector()
        self.process.buildActions.percentile_95.vectorKey = "flatDistribution_95.0"

        self.process.buildActions.percentile_100 = LoadVector()
        self.process.buildActions.percentile_100.vectorKey = "flatDistribution_100.0"

        # self.process.calculateActions.mag50 = Mag50Action()
        # self.process.calculateActions.mag50.vectorKey = "{band}_mag_ref"
        # self.process.calculateActions.mag50.matchDistanceKey = "matchDistance"

        self.produce.plot = PercentilePlot()
        # self.produce.metric.units = {"mag50": "mag"}
        # self.produce.metric.newNames = {"mag50": "{band}_mag50"}
