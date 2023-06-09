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
from __future__ import annotations

import os
from unittest import TestCase, main

import lsst.utils.tests
from lsst.analysis.tools.atools import SkyObjectHistPlot
from lsst.analysis.tools.interfaces import AnalysisBaseConfig
from lsst.analysis.tools.tasks.objectTableTractAnalysis import ObjectTableTractAnalysisConfig
from lsst.analysis.tools.tasks.reconstructor import getPlotDatasetTypeNames
from lsst.daf.butler import Butler, DatasetType
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.pex.config import Config

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class TestGetPlotDatasetTypeNames(TestCase):
    """Test plot dataset type names are correctly queried by the butler."""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)
        Butler.makeRepo(self.root)
        self.butler = Butler(self.root, run="test_run")
        # No dimensions in dataset type so we don't have to worry about
        # inserting dimension data or defining data IDs.
        self.datasetTypePlot = DatasetType(
            "plot1Label_config",
            dimensions=(),
            storageClass="Config",
            universe=self.butler.dimensions,
        )
        self.butler.registry.registerDatasetType(self.datasetTypePlot)
        self.datasetTypePlot = DatasetType(
            "plot2Label_config",
            dimensions=(),
            storageClass="Config",
            universe=self.butler.dimensions,
        )
        self.butler.registry.registerDatasetType(self.datasetTypePlot)
        self.datasetTypeBase = DatasetType(
            "baseLabel_config",
            dimensions=(),
            storageClass="Config",
            universe=self.butler.dimensions,
        )
        self.butler.registry.registerDatasetType(self.datasetTypeBase)
        self.datasetTypeOther = DatasetType(
            "otherLabel_config",
            dimensions=(),
            storageClass="Config",
            universe=self.butler.dimensions,
        )
        self.butler.registry.registerDatasetType(self.datasetTypeOther)

    def tearDown(self):
        removeTestTempDir(self.root)

    def testGetPlotNames(self):
        """Test that we can get plot dataset type names from a given collection
        in a butler. Two queries are made: one across the entire collection,
        and one for only dataset type names prefixed with a certain label.
        """
        plot1Config = ObjectTableTractAnalysisConfig()
        plot1Config.plots.skyObjectHistPlot1 = SkyObjectHistPlot()
        plot2Config = ObjectTableTractAnalysisConfig()
        plot2Config.plots.skyObjectHistPlot2 = SkyObjectHistPlot()
        baseConfig = AnalysisBaseConfig()
        baseConfig.connections.outputName = "baseName"
        otherConfig = Config()
        _ = self.butler.put(plot1Config, "plot1Label_config")
        _ = self.butler.put(plot2Config, "plot2Label_config")
        _ = self.butler.put(baseConfig, "baseLabel_config")
        _ = self.butler.put(otherConfig, "otherLabel_config")
        plots_all = getPlotDatasetTypeNames(self.butler, "test_run")
        plots_label = getPlotDatasetTypeNames(self.butler, "test_run", "plot1Label")
        for plot in plots_all:
            self.assertIn("skyObjectHistPlot", plot)
        for plot in plots_label:
            self.assertIn("skyObjectHistPlot1", plot)


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
