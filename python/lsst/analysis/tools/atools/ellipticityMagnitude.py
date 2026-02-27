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

from ..actions.vector import LoadVector
from ..actions.vector.mathActions import ConstantValue, DivideVector, SubtractVector
from .actionMagnitudeScatterPlot import ActionMagnitudeScatterPlot

__all__ = (
    "EllipticityMagnitudePlotBase",
    "ExponentialEllipticityMagnitudePlot",
    "SersicEllipticityMagnitudePlot",
)


class EllipticityMagnitudePlotBase(ActionMagnitudeScatterPlot):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = (0.0, 1.0)
        self.produce.plot.legendLocation = "upper right"


class ExponentialEllipticityMagnitudePlot(EllipticityMagnitudePlotBase):
    def setDefaults(self):
        super().setDefaults()
        self.action_column = SubtractVector(
            actionA=ConstantValue(value=1),
            actionB=DivideVector(
                actionA=LoadVector(vectorKey="exponential_reff_minor"),
                actionB=LoadVector(vectorKey="exponential_reff_major"),
            ),
        )
        self.key_y = "axrat_exponential"
        self.produce.plot.yAxisLabel = "Exponential ellipticity"


class SersicEllipticityMagnitudePlot(EllipticityMagnitudePlotBase):
    def setDefaults(self):
        super().setDefaults()
        self.action_column = SubtractVector(
            actionA=ConstantValue(value=1),
            actionB=DivideVector(
                actionA=LoadVector(vectorKey="sersic_reff_minor"),
                actionB=LoadVector(vectorKey="sersic_reff_major"),
            ),
        )
        self.key_y = "axrat_sersic"
        self.produce.plot.yAxisLabel = "Sérsic ellipticity"
        self.produce.plot.yLims = (0.0, 1.0)
