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
from lsst.analysis.tools.actions.vector import LoadVector
from lsst.analysis.tools.interfaces import AnalysisTool, JointAction, PlotAction

# BaseMetricAction is used to produce a real Measurement without needing to
# hand-roll a MetricAction subclass.
from lsst.analysis.tools.interfaces._stages import BaseMetricAction
from lsst.verify import Measurement


class NullPlot(PlotAction):
    def __call__(self, data, **kwargs):
        return plt.figure()


class ProducePlotsToolTestCase(TestCase):
    """Test the `AnalysisTool.doProducePlots` switch."""

    def _makeJointTool(self) -> AnalysisTool:
        tool = AnalysisTool()
        tool.prep.keysToLoad = ["value"]
        tool.process.buildActions.value = LoadVector(vectorKey="value")
        joint = JointAction(metric=BaseMetricAction(), plot=NullPlot())
        joint.metric.units = {"value": "ct"}
        tool.produce = joint
        return tool

    def testJointActionSkipsPlotOnly(self):
        """With doProducePlots=False, the metric is still produced but the
        plot is not.
        """
        tool = self._makeJointTool()
        data = {"value": 5}

        resultsWithPlots = tool(data)
        self.assertTrue(any(isinstance(v, Measurement) for v in resultsWithPlots.values()))
        self.assertTrue(any(hasattr(v, "savefig") for v in resultsWithPlots.values()))

        tool.doProducePlots = False
        resultsNoPlots = tool(data)
        self.assertTrue(any(isinstance(v, Measurement) for v in resultsNoPlots.values()))
        self.assertFalse(any(hasattr(v, "savefig") for v in resultsNoPlots.values()))

    def testPurePlotActionProducesNothing(self):
        """An AnalysisTool whose produce action is only a PlotAction (no
        metric) produces no results at all when doProducePlots is False.
        """
        tool = AnalysisTool()
        tool.produce = NullPlot()

        self.assertTrue(tool({}))

        tool.doProducePlots = False
        self.assertEqual(tool({}), {})

    def testGetOutputNamesEmptyWhenDisabled(self):
        tool = self._makeJointTool()
        self.assertNotEqual(tuple(tool.getOutputNames()), tuple())

        tool.doProducePlots = False
        self.assertEqual(tuple(tool.getOutputNames()), tuple())


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
