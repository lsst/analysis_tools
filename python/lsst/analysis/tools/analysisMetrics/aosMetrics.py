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

__all__ = ("ZernikesFwhmMetricLsst", "ZernikesFwhmMetricLatiss")

from ..actions.vector import CalcFwhmZernikesLsst, CalcFwhmZernikesLatiss
from ..interfaces import AnalysisMetric


class ZernikesFwhmMetricLsst(AnalysisMetric):
    """Calculate FWHM from Zernikes for LsstCam."""

    def setDefaults(self):
        super().setDefaults()

        # Calculate FWHM from Zernike Coefficients
        self.process.calculateActions.ZernikesFwhmMetric = CalcFwhmZernikesLsst(
            vectorKey="zernikeEstimateAvg"
        )

        self.produce.units = {"ZernikesFwhmMetric": "arcsec"}


class ZernikesFwhmMetricLatiss(AnalysisMetric):
    """Calculate FWHM from Zernikes for LAtiss."""

    def setDefaults(self):
        super().setDefaults()

        # Calculate FWHM from Zernike Coefficients
        self.process.calculateActions.ZernikesFwhmMetric = CalcFwhmZernikesLatiss(
            vectorKey="zernikeEstimateAvg"
        )

        self.produce.units = {"ZernikesFwhmMetric": "arcsec"}
