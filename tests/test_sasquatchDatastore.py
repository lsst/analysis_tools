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

import datetime
import os
import unittest
from unittest.mock import patch

import astropy.units as u

import lsst.daf.butler.tests as butlerTests
from lsst.analysis.tools.interfaces import MetricMeasurementBundle
from lsst.analysis.tools.interfaces.datastore import (
    DispatchableMetricBundle,
    MetricMeasurement,
    MetricName,
    SasquatchDispatcher,
)
from lsst.daf.butler import CollectionType, Config, DatasetTypeNotSupportedError
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.verify import Measurement

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class DuckMetricBundle:
    """A metric bundle from outside this package.

    Implements the dispatcher interface without inheriting from
    `MetricMeasurementBundle`, standing in for the type a site other
    than Rubin might configure this datastore to accept.
    """

    def __init__(self, data=None):
        self._data = dict(data or {})
        self.metricNamePrefix = ""

    def items(self):
        return self._data.items()


class SasquatchDatastoreTest(unittest.TestCase):
    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

        config = Config()
        config["datastore", "cls"] = "lsst.analysis.tools.interfaces.datastore.SasquatchDatastore"
        config["datastore", "restProxyUrl"] = "https://example.com/sasquatch-rest-proxy"
        config["storageClasses", "DuckMetricBundle"] = {"pytype": f"{__name__}.DuckMetricBundle"}

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

    def test_put_duck_typed_bundle(self):
        """A type meeting the protocol dispatches whatever its storage
        class is called, so sites outside Rubin can configure their own.
        """
        butlerTests.addDatasetType(
            self.butler, "DuckMetrics", {"instrument", "visit", "detector"}, "DuckMetricBundle"
        )
        m = Measurement("nopackage.fancyMetric", 42.2 * u.s)
        bundle = DuckMetricBundle({"m": [m]})

        with patch.object(SasquatchDispatcher, "dispatchRef") as mock_method:
            self.butler.put(bundle, "DuckMetrics", run="run1", instrument="DummyCam", visit=42, detector=2)

        mock_method.assert_called()
        self.assertIs(mock_method.call_args[0][0], bundle)

    def test_rejects_non_conforming_measurements(self):
        """A bundle holding things that are not measurements is rejected
        before dispatch rather than failing inside the dispatcher.
        """
        butlerTests.addDatasetType(
            self.butler, "DuckMetrics", {"instrument", "visit", "detector"}, "DuckMetricBundle"
        )
        bundle = DuckMetricBundle({"m": ["not a measurement"]})

        with patch.object(SasquatchDispatcher, "dispatchRef") as mock_method:
            with self.assertRaises(DatasetTypeNotSupportedError):
                self.butler.put(
                    bundle, "DuckMetrics", run="run1", instrument="DummyCam", visit=42, detector=2
                )

        mock_method.assert_not_called()

    def test_put_bundle_with_no_measurements(self):
        """A bundle with nothing to measure has nothing to check and is
        dispatched unchanged.
        """
        butlerTests.addDatasetType(
            self.butler, "DuckMetrics", {"instrument", "visit", "detector"}, "DuckMetricBundle"
        )
        bundle = DuckMetricBundle({"empty": []})

        with patch.object(SasquatchDispatcher, "dispatchRef") as mock_method:
            self.butler.put(bundle, "DuckMetrics", run="run1", instrument="DummyCam", visit=42, detector=2)

        mock_method.assert_called()

    def test_rejects_unsupported_dataset_type(self):
        """A storage class the configuration accepts is still rejected when
        its Python type does not provide the dispatcher's interface.
        """
        butlerTests.addDatasetType(
            self.butler, "NotAMetric", {"instrument", "visit", "detector"}, "StructuredDataDict"
        )
        with self.assertRaisesRegex(DatasetTypeNotSupportedError, "metricNamePrefix"):
            self.butler.put({"a": 1}, "NotAMetric", run="run1", instrument="DummyCam", visit=42, detector=2)

    def test_explicit_timestamp_version(self):
        dispatcher = SasquatchDispatcher("http://test.local", "na")
        bundle = MetricMeasurementBundle()
        bundle.timestamp_version = "explicit_timestamp"
        # verify this raises with no specified time
        with self.assertRaises(ValueError):
            dispatcher._handleTimes({}, bundle, "localRun")
        # verify this raise with a date that can't be parsed
        bundle.timestamp_version = "explicit_timestamp:123233"
        with self.assertRaises(ValueError):
            dispatcher._handleTimes({}, bundle, "localRun")
        # verify that a correct time gets parsed
        bundle.timestamp_version = "explicit_timestamp:20230728T165102Z"
        meta = {}
        dispatcher._handleTimes(meta, bundle, "localRun")
        dt = datetime.datetime(2023, 7, 28, 16, 51, 2, tzinfo=datetime.UTC)
        self.assertEqual(meta["timestamp"], dt.timestamp())


class DispatchProtocolTest(unittest.TestCase):
    """Tests for the structural types describing the dispatcher's interface."""

    def test_bundle_accepts_metric_measurement_bundle(self):
        self.assertIsInstance(MetricMeasurementBundle(), DispatchableMetricBundle)

    def test_bundle_rejects_plain_dict(self):
        """A `dict` has ``items`` but no ``metricNamePrefix``."""
        self.assertNotIsInstance({"a": 1}, DispatchableMetricBundle)

    def test_measurement_accepts_verify_measurement(self):
        measurement = Measurement("nopackage.fancyMetric", 42.2 * u.s)
        self.assertIsInstance(measurement, MetricMeasurement)
        self.assertIsInstance(measurement.metric_name, MetricName)

    def test_measurement_rejects_unrelated_object(self):
        self.assertNotIsInstance(object(), MetricMeasurement)


if __name__ == "__main__":
    unittest.main()
