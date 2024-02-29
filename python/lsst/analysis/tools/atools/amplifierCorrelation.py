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
    "ImageCorrelationPlot",
    "OverscanCorrelationPlot",
)

from ..actions.plot.matrixPlot import MatrixPlot
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool


class BaseCorrelationPlot(AnalysisTool):
    """Base class for correlation plots of amplifier biases."""

    # Do not iterate over multiple bands in a parameterized manner because this
    # AnalysisTool does not support band as a name parameter.
    parameterizedBand: bool = False

    def setupPlot(self, matrixKey: str, title: str):
        """
        Set up a matrix plot action.

        Parameters
        ----------
        matrixKey : `str`
            Key for the matrix being plotted.
        title : `str`
            Title for the plot.

        Returns
        -------
        action : `~lsst.analysis.tools.actions.plot.matrixPlot.MatrixPlot`
            The resulting plot action.
        """
        action = MatrixPlot()
        action.matrixKey = matrixKey
        action.component1Key = "ampComp"
        action.component2Key = "ampName"
        action.componentGroup1Key = "detectorComp"
        action.componentGroup2Key = "detector"
        action.setPositionsAtPixelBoundaries = True
        action.hideMinorTicks = ["x", "y"]
        action.tickLabelsFontSize = 6
        action.title = title
        action.titleFontSize = 9
        action.xAxisLabel = "Amplifiers in detector "
        action.yAxisLabel = "Amplifiers in detector "
        action.colorbarLabel = "Correlation value"
        action.colorbarLabelFontSize = 9
        action.colorbarTickLabelFontSize = 7.5
        return action

    def setDefaults(self):
        super().setDefaults()

        # Initialize plot with specific parameters.
        action = self.setupPlot(matrixKey=self.matrixKey, title=self.plotTitle)

        # Load the relevant columns.
        for key in (
            action.matrixKey,
            action.component1Key,
            action.component2Key,
            action.componentGroup1Key,
            action.componentGroup2Key,
        ):
            setattr(self.process.buildActions, key, LoadVector(vectorKey=key))

        # Provide the plot action to the AnalysisTool.
        self.produce.plot = action


class ImageCorrelationPlot(BaseCorrelationPlot):
    """Plot image correlation of amplifier biases."""

    matrixKey = "imageCorr"
    plotTitle = "Image correlations"


class OverscanCorrelationPlot(BaseCorrelationPlot):
    """Plot overscan correlation of amplifier biases."""

    matrixKey = "overscanCorr"
    plotTitle = "Overscan correlations"
