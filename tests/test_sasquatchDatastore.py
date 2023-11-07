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

import os
import unittest
from unittest.mock import patch

import astropy.units as u
import lsst.daf.butler.tests as butlerTests
from lsst.analysis.tools.interfaces import MetricMeasurementBundle
from lsst.analysis.tools.interfaces.datastore import SasquatchDispatcher
from lsst.daf.butler import CollectionType, Config
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.verify import Measurement

TESTDIR = os.path.abspath(os.path.dirname(__file__))
CONFIG_FILE = os.path.join(TESTDIR, "config", "butler-sasquatch.yaml")


class SasquatchDatastoreTest(unittest.TestCase):
    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

        config = Config()
        config["datastore", "cls"] = "lsst.analysis.tools.interfaces.datastore.SasquatchDatastore"
        config["datastore", "restProxyUrl"] = "https://example.com/sasquatch-rest-proxy"

        dataIds = {
            "instrument": ["DummyCam"],
            "physical_filter": ["d-r"],
            "visit": [42, 43, 44],
            "detector": [1, 2, 3],
        }
        self.butler = butlerTests.makeTestRepo(self.root, dataIds, config=config)

        butlerTests.addDatasetType(
            self.butler, "Metrics", {"instrument", "visit", "detector"}, "MetricMeasurementBundle"
        )
        self.butler.registry.registerCollection("run1", CollectionType.RUN)

    def tearDown(self):
        removeTestTempDir(self.root)

    def test_put(self):
        """Simple test for put method."""
        m = Measurement("nopackage.fancyMetric", 42.2 * u.s)
        bundle = MetricMeasurementBundle({"m": [m]})

        # Patch dispatcher method to check parameters.
        with patch.object(SasquatchDispatcher, "dispatchRef") as mock_method:
            self.butler.put(bundle, "Metrics", run="run1", instrument="DummyCam", visit=42, detector=2)

        mock_method.assert_called()
        self.assertIs(mock_method.call_args[0][0], bundle)


if __name__ == "__main__":
    unittest.main()
