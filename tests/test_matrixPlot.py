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

import unittest

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import lsst.utils.tests
from lsst.analysis.tools.actions.plot.matrixPlot import GuideLinesConfig, MatrixPlot
from lsst.analysis.tools.actions.vector.vectorActions import LoadVector
from lsst.analysis.tools.interfaces import AnalysisTool

# No display needed.
matplotlib.use("Agg")


class MatrixPlotTaskTestCase(lsst.utils.tests.TestCase):
    """Test to catch basic errors or logical inconsistencies in matrix plot
    generation.
    """

    def setUp(self):
        # Configure the plot action.
        action = MatrixPlot()
        action.matrixKey = "custom_key"
        action.matrixOrigin = "lower"
        action.xAxisLabel = "X axis label"
        action.yAxisLabel = "Y axis label"
        action.colorbarLabel = "Colorbar label"
        action.title = "Title"
        action.titleFontSize = 7
        action.figsize = [4, 4]
        action.dpi = 310
        action.guideLines["x"] = GuideLinesConfig()
        action.guideLines["x"].lines = {2.5: "x=2.5", 6: ""}
        action.guideLines["x"].color = "red"
        action.guideLines["x"].linestyle = "--"
        action.guideLines["y"] = GuideLinesConfig()
        action.guideLines["y"].lines = {3.5: "Another label @ y = 3.5"}
        action.guideLines["y"].color = "white"
        action.guideLines["y"].outlineColor = "orange"
        action.guideLines["y"].linestyle = ":"
        action.xAxisTickValues = [0, 2, 4, 6, 8, 10]
        action.xAxisTickLabels = {1: "X0", 3: "X1", 5: "X2", 7: "X3", 9: "X4"}
        action.yAxisTickValues = [0, 3, 6, 8, 10]
        action.yAxisTickLabels = {1.5: "Y0", 4.5: "Y1", 7: "Y2", 9: "Y3"}
        action.setPositionsAtPixelBoundaries = True
        action.hideMinorTicks = ["y"]

        # Set up `AnalysisTool` for the plot action.
        self.plot = AnalysisTool()
        self.plot.produce.plot = action

        # Mock up the data.
        np.random.seed(1905)
        matrix = np.random.rand(10, 10)
        self.data = {"custom_key": matrix}

        # Load the relevant column and finalize the plot.
        self.plot.process.buildActions.custom_key = LoadVector(vectorKey="custom_key")
        self.plot.finalize()

        # Set up the plotInfo dictionary.
        self.plotInfo = {key: "test" for key in ("plotName", "run", "tableName")}
        self.plotInfo["bands"] = []

    def tearDown(self):
        del self.plot
        del self.data
        del self.plotInfo

    def test_MatrixPlotTask(self):
        plt.rcParams.update(plt.rcParamsDefault)

        # Run the `AnalysisTool` for the plot action and get the result.
        result = self.plot(data=self.data, skymap=None, plotInfo=self.plotInfo)

        # Get the type of the plot.
        plot_type = type(self.plot.produce.plot).__name__

        # Assert that the `plot_type` key exists in the result dictionary.
        self.assertTrue(plot_type in result, msg=f"{plot_type} key is missing in the result dictionary")

        # Unpack the figure object from the dictionary.
        fig = result[plot_type]

        # Check that the returned object is indeed a matplotlib figure.
        self.assertTrue(isinstance(fig, plt.Figure), msg="The output is not a matplotlib figure.")

        # Assert that the figure has the correct size and dpi.
        figsize = self.plot.produce.plot.figsize
        dpi = self.plot.produce.plot.dpi
        self.assertTrue(
            all(fig.get_size_inches() == figsize),
            msg=f"Figure size is not {figsize}, it is {fig.get_size_inches()}",
        )
        self.assertTrue(fig.dpi == dpi, msg=f"Figure dpi is not {dpi}, it is {fig.dpi}")

        # Assert that the figure has 2 axes.
        self.assertTrue(len(fig.axes) == 2, f"Figure does not have 2 axes, it has {len(fig.axes)}")


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
