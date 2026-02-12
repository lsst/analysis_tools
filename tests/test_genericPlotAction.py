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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from unittest import TestCase, main

import matplotlib.pyplot as plt

import lsst.utils.tests
from lsst.analysis.tools.atools.genericPlotAction import StructPlotAction
from lsst.analysis.tools.interfaces import AnalysisTool, JointAction, NoMetric, PlotAction


class NullPlot(PlotAction):
    def __call__(self, data, **kwargs):
        fig = plt.figure()
        return fig


class CustomPlot(NullPlot):
    def getPlotType(self) -> str:
        return "Custom"


class StructPlotActionTestCase(TestCase):
    """Test the generic StructPlotAction"""

    def setUp(self) -> None:
        super().setUp()
        action_plot = StructPlotAction()
        action_plot.actions.hist = NullPlot()
        action_plot.actions.xy = CustomPlot()
        action_joint = JointAction(metric=NoMetric(), plot=action_plot)
        self.action_joint = action_joint
        tool = AnalysisTool()
        tool.produce = action_joint
        self.tool = tool

    def testTool(self) -> None:
        results = self.tool({})
        assert tuple(results.keys()) == ("hist_NullPlot", "xy_Custom")


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
