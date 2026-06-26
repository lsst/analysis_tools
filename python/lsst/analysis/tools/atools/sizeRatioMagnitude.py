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

from lsst.pex.config import Field

from ..actions.vector import LoadVector
from ..actions.vector.mathActions import ConstantValue, DivideVector, MultiplyVector
from .actionMagnitudeScatterPlot import ActionMagnitudeScatterPlot

__all__ = ("SizeRatioMagnitudePlot",)


class SizeRatioMagnitudePlot(ActionMagnitudeScatterPlot):
    factor_mult = Field[float](doc="Factor to multiply size ratio by", default=1.0)
    key_divisor = Field[str](doc="Key for the size field to use as the divisor")

    def setDefaults(self):
        super().setDefaults()
        self.action_column = DivideVector(
            actionA=MultiplyVector(
                actionA=ConstantValue(value=1.0),
                actionB=LoadVector(vectorKey="sersic_reff_major"),
            ),
            actionB=LoadVector(vectorKey="placeholder"),
        )
        self.key_y = "size_ratio"
        self.produce.plot.yAxisLabel = "Size ratio"
        self.produce.plot.yLims = (0.0, 1.0)

    def finalize(self):
        super().finalize()
        self.action_column.actionA.actionA.value = self.factor_mult
        self.action_column.actionB.vectorKey = self.key_divisor
