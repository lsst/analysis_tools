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
    "DiaDipoleOrientationVsParallacticAngle",
    "DiaDipoleOrientationVsDiaKernelShift",
    "DiaSourcesGoodVsBadRatioMetric",
    "NumDiaSourcesSelectionMetric",
    "NumDipolesMetric",
    "NumGoodDiaSourcesMetrics",
)

from lsst.pex.config import Field

from ..actions.scalar import CountAction, DivideScalar
from ..actions.vector import DownselectVector, FlagSelector, GoodDiaSourceSelector, LoadVector
from ..interfaces import AnalysisTool


class NumGoodDiaSourcesMetrics(AnalysisTool):
    """Calculate the number of DIA Sources that do not have known
    bad/quality flags set to true, and also calculate the ratio of
    counts of non-flagged sources to all sources.
    """

    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        # filter for and count the number of dia sources that don't have flags
        self.process.filterActions.goodDiaSources = DownselectVector(
            vectorKey="parentDiaSourceId", selector=GoodDiaSourceSelector()
        )
        self.process.calculateActions.numGoodDiaSources = CountAction(vectorKey="goodDiaSources")

        # Count the total number of dia sources:
        self.process.calculateActions.numAllDiaSources = CountAction(vectorKey="parentDiaSourceId")

        # And calculate the ratio of good-to-all counts
        self.process.calculateActions.ratioGoodToAllDiaSources = DivideScalar(
            actionA=self.process.calculateActions.numGoodDiaSources,
            actionB=self.process.calculateActions.numAllDiaSources,
        )

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {
            "numAllDiaSources": "ct",
            "numGoodDiaSources": "ct",
            "ratioGoodToAllDiaSources": "",
        }


class NumDipolesMetric(AnalysisTool):
    """Calculate the number of dipoles with NaN values excluded."""

    def setDefaults(self):
        super().setDefaults()

        # select all diaSources flagged as dipole
        self.prep.selectors.flags = FlagSelector(selectWhenTrue=["isDipole"])

        # count the number of dipoles
        self.process.buildActions.numDipoles = CountAction(vectorKey="isDipole")

        # the units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"numDipoles": "ct"}


class NumDiaSourcesSelectionMetric(AnalysisTool):
    """Count the number of DIA Sources for a given threshold."""

    metricName = Field[str](doc="Name to use for output metric")

    def setDefaults(self):
        super().setDefaults()

        # Count dia sources with reliability lower than the threshold
        self.process.calculateActions.countingAction = CountAction

        # The units for the quantity (count, an astropy quantity)
        self.produce.metric.units = {"countingAction": "ct"}

    def finalize(self):
        self.produce.metric.newNames = {"countingAction": self.metricName}


class DiaSourcesGoodVsBadRatioMetric(AnalysisTool):
    """Calculate the ratio of 'good' vs 'bad' DIA Sources."""

    def setDefaults(self):
        super().setDefaults()

        # Count dia sources with reliability higher than the threshold
        self.process.buildActions.numDiaSourcesHighReliability = CountAction(
            op="gt", threshold=0.9, vectorKey="reliability"
        )

        # Count dia sources with reliability lower than the threshold
        self.process.buildActions.numDiaSourcesLowReliability = CountAction(
            op="lt", threshold=0.1, vectorKey="reliability"
        )

        # Calculate ratio of good vs bad DIA Sources
        self.process.calculateActions.DiaSourcesGoodVsBadRatio = DivideScalar(
            actionA=self.process.buildActions.numDiaSourcesHighReliability,
            actionB=self.process.buildActions.numDiaSourcesLowReliability,
        )

        # The units for the quantity (dimensionless, an astropy quantity)
        self.produce.metric.units = {"DiaSourcesGoodVsBadRatio": ""}


class DiaSourcesDipoleOrientationVsParallacticAngleMetric(AnalysisTool):
    """Calculate the dipole orientation w.r.t. parallactic angle.
    
    Calculates the von Mises mu, kappa

    Defines the reported quantities as
    mean = mu
    sigma = sqrt(1 / kappa)
    """

    def setDefaults(self):
        super().setDefaults()

        # What sets the capitalization convention here?
        self.process.buildActions.parallacticAngle = CalculateParallacticAngle(
            diff
        )
        # this is the orientation in detector, y from x coordinates.
        self.process.buildActions.dipoleOrientation = LoadVector(vectorKey="ip_diffim_DipoleFit_orientation")
        # Are we allowed to do arithmetic operations directly?
        self.process.buildActions.dipoleOrientationVsParallacticAngle = 
            self.process.buildActions.dipoleOrientation - \
            self.process.buildActions.parallacticAngle
        self.process.calculateActions.DiaSourcesDipoleOrientationVsParallacticAngle = fitVonMises(
            dipole_orientation_vs_parallactic_angle
        )

        # DIA Kernel
        self.process.buildActions.diaKernelShift = CalculateDiaKernelShift(
            diaKernel
        )
        self.process.buildActions.dipoleOrientationVsDiaShiftAngle = self.process.buildActions.dipoleOrientaiton - self.process.buildActions.diaKernelShift
        self.process.buildActions.dipoleOrientationVsDiaShift = fitVonMises(
            dipole_orientation_vs_dia_shift_angle
        )
        self.process.calculateActions.DiaSourcesDipoleOrientationVsParallacticAngle = fitVonMises(
            dipole_orientation_vs_parallactic_angle,
        )
        self.produce.metric.units = {
            "diaSourcesDipoleOrientationVsParallacticAngleMean": "rad",
            "diaSourcesDipoleOrientationVsParallacticAngleSigma": "rad",
            "diaSourcesDipoleOrientationVsDiaKernelShiftAngleMean": "rad",
            "diaSourcesDipoleOrientationVsDiaKernelShiftAngleSigma": "rad",
        }
