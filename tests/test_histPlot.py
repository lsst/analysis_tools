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

import lsst.utils.tests
from lsst.analysis.tools.actions.plot.histPlot import HistPanel, HistPlot

matplotlib.use("Agg")


class HistPlotLayoutTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.plot = HistPlot()
        for panel_name in ("panel_a", "panel_b", "panel_c"):
            self.plot.panels[panel_name] = HistPanel()
            self.plot.panels[panel_name].hists = {f"{panel_name}_hist": panel_name}

    def test_default_layout_uses_two_columns(self):
        fig = plt.figure()

        axes, ncols, nrows = self.plot._makeAxes(fig)

        self.assertEqual(ncols, 2)
        self.assertEqual(nrows, 2)
        self.assertEqual(len(axes), 3)

    def test_panels_per_row_one_stacks_vertically(self):
        self.plot.panelsPerRow = 1
        fig = plt.figure()

        axes, ncols, nrows = self.plot._makeAxes(fig)

        self.assertEqual(ncols, 1)
        self.assertEqual(nrows, 3)
        self.assertEqual(len(axes), 3)
        self.assertTrue(all(ax.get_position().x0 == axes[0].get_position().x0 for ax in axes[1:]))
        self.assertGreater(axes[0].get_position().y0, axes[1].get_position().y0)
        self.assertGreater(axes[1].get_position().y0, axes[2].get_position().y0)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
