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

import lsst.utils.tests
import matplotlib.pyplot as plt
from astropy.table import Table
from lsst.analysis.tools.atools import diaSourceTableTractMetrics as diaAtool


class DiaSourceTableTractTest(lsst.utils.tests.TestCase):
    """Test to see if DiaSources are counted as expected
    and a plot is generated without falling over."""

    def setUp(self):
        testFile = "tests/diaSourceTable_tract_test.ecsv"
        self.data = Table.read(testFile)

    def test_metrics(self):
        """Test that metrics have the expected values from the test data."""
        NumDiaSources = diaAtool.NumDiaSourcesMetric()
        NumDiaSources.finalize()
        goodDiaSourceCount = NumDiaSources(self.data)
        self.assertEqual(goodDiaSourceCount["numDiaSources"].quantity.value, 344)

        NumStreakDiaSources = diaAtool.NumStreakDiaSourcesMetric()
        NumStreakDiaSources.finalize()
        streakDiaSourceCount = NumStreakDiaSources(self.data)
        self.assertEqual(streakDiaSourceCount["numStreakDiaSources"].quantity.value, 3)

        NumStreakCenterDiaSources = diaAtool.NumStreakCenterDiaSourcesMetric()
        NumStreakCenterDiaSources.finalize()
        streakCenterDiaSourceCount = NumStreakCenterDiaSources(self.data)
        self.assertEqual(streakCenterDiaSourceCount["numStreakCenterDiaSources"].quantity.value, 1)

    def test_plot(self):
        """Test that a plot is created from the test data."""
        plot = diaAtool.PlotStreakDiaSources()
        plot.finalize()
        result = plot(self.data)
        self.assertTrue(isinstance(result["DiaSkyPlot"], plt.Figure))


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
