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

"""Tests for the calibQuantityProfile tools on cp_verify results tables.

These catch rows with no panel key (the per-detector and per-exposure rows,
masked in the amplifier column) breaking the repack, and with it
analyzeFlatDetCore's flatTestsByDate plot.
"""

import unittest

import matplotlib
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from matplotlib.figure import Figure

import lsst.utils.tests
from lsst.analysis.tools.atools.calibQuantityProfile import CalibAmpScatterTool, SingleValueRepacker

matplotlib.use("Agg")

AMP_NAMES = [f"C{i}{j}" for i in (0, 1) for j in range(8)]
DETECTORS = (90, 91)
MJD = 61224.4


def makeResultsTable() -> Table:
    """Make a table shaped like cp_verify's ``verifyFlatResults``.

    Returns
    -------
    table : `astropy.table.Table`
        One row per amplifier per detector, plus a per-detector row for each
        detector and a per-exposure row. Those last have no amplifier, so are
        masked in the ``amplifier`` column.
    """
    nAmpRows = len(DETECTORS) * len(AMP_NAMES)
    ampRows = Table(
        {
            "detector": [detector for detector in DETECTORS for _ in AMP_NAMES],
            "amplifier": AMP_NAMES * len(DETECTORS),
            "mjd": [MJD] * nAmpRows,
            "FLAT_VERIFY_NOISE": [True] * nAmpRows,
        }
    )
    detectorRows = Table(
        {"detector": list(DETECTORS), "mjd": [MJD] * len(DETECTORS), "FLAT_DET_SCATTER": [0.01, 0.01]}
    )
    exposureRow = Table({"detector": [DETECTORS[0]], "flat_SCATTER": [0.001]})
    return vstack([ampRows, detectorRows, exposureRow])


class SingleValueRepackerTestCase(lsst.utils.tests.TestCase):
    def testSkipsRowsWithoutPanelKey(self) -> None:
        # The rows with no amplifier used to raise
        # "TypeError: unhashable type: 'MaskedConstant'", so
        # analyzeFlatDetCore never wrote its plot.
        table = makeResultsTable()
        # Without masked rows in the fixture this test pins nothing.
        self.assertTrue(table["amplifier"].mask.any())
        repacker = SingleValueRepacker(panelKey="amplifier", dataKey="mjd", quantityKey="FLAT_VERIFY_NOISE")
        repacked = repacker(table)
        self.assertEqual(set(repacked), {f"{amp}{suffix}" for amp in AMP_NAMES for suffix in ("", "_x")})
        for amp in AMP_NAMES:
            self.assertEqual(repacked[f"{amp}_x"], [MJD] * len(DETECTORS))
            self.assertEqual(repacked[amp], [True] * len(DETECTORS))


class CalibAmpScatterToolTestCase(lsst.utils.tests.TestCase):
    def testPlotsResultsTable(self) -> None:
        # The whole tool, configured as analyzeFlatDetCore's flatTestsByDate
        # in cpCore.yaml: a masked amplifier anywhere along the way must not
        # stop the plot being made.
        tool = CalibAmpScatterTool()
        tool.prep.quantityKey = "FLAT_VERIFY_NOISE"
        tool.finalize()
        results = tool(makeResultsTable())
        self.assertIsInstance(results["GridPlot"], Figure)
        plt.close(results["GridPlot"])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
