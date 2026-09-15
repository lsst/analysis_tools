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
from unittest.mock import MagicMock, patch

import astropy.units as u
import requests

import lsst.daf.butler.tests as butlerTests
from lsst.analysis.tools.interfaces import MetricMeasurementBundle
from lsst.analysis.tools.interfaces.datastore import SasquatchDispatcher, SasquatchDispatchFailure
from lsst.analysis.tools.interfaces.datastore._dispatcher import DEFAULT_TIMEOUT
from lsst.daf.butler import CollectionType, Config
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.verify import Measurement

TESTDIR = os.path.abspath(os.path.dirname(__file__))
CONFIG_FILE = os.path.join(TESTDIR, "config", "butler-sasquatch.yaml")

# Where the dispatcher looks up its HTTP session factory, for patching.
HTTP_CLIENT = "lsst.analysis.tools.interfaces.datastore._dispatcher.http_client"
DATASTORE_LOGGER = "lsst.analysis.tools.interfaces.datastore._sasquatchDatastore"


def _fakeResponse(status: int = 200, json=None, jsonError: Exception | None = None) -> MagicMock:
    """Make a stand-in for a `requests.Response`.

    ``raise_for_status`` raises for 4xx/5xx statuses, and ``json()`` either
    returns ``json`` or raises ``jsonError``.
    """
    response = MagicMock()
    response.status_code = status
    if status >= 400:
        response.raise_for_status.side_effect = requests.HTTPError(f"HTTP {status}", response=response)
    if jsonError is not None:
        response.json.side_effect = jsonError
    else:
        response.json.return_value = json
    return response


def _fakeSession(getResponse: MagicMock | None = None, postResponse: MagicMock | None = None) -> MagicMock:
    """Make a stand-in for a `requests.Session` returning canned responses."""
    session = MagicMock()
    session.get.return_value = getResponse if getResponse is not None else _fakeResponse()
    session.post.return_value = postResponse if postResponse is not None else _fakeResponse()
    return session


class SasquatchDatastoreTest(unittest.TestCase):
    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

        config = Config()
        config["datastore", "cls"] = "lsst.analysis.tools.interfaces.datastore.SasquatchDatastore"
        config["datastore", "restProxyUrl"] = "https://example.com/sasquatch-rest-proxy"
        config["datastore", "timeout"] = 5.5

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

    def test_put_survives_unexpected_dispatch_error(self):
        """An unexpected exception from the dispatcher must not fail the put.

        Publishing to Sasquatch is best-effort, so a broken proxy (or a bug
        in preparing the records) must never take down the pipeline task
        that wrote the metric bundle. The failure is logged instead.
        """
        m = Measurement("nopackage.fancyMetric", 42.2 * u.s)
        bundle = MetricMeasurementBundle({"m": [m]})

        with patch.object(SasquatchDispatcher, "dispatchRef", side_effect=KeyError("data")):
            with self.assertLogs(DATASTORE_LOGGER, level="ERROR") as logs:
                ref = self.butler.put(
                    bundle, "Metrics", run="run1", instrument="DummyCam", visit=43, detector=2
                )

        self.assertIsNotNone(ref)
        self.assertTrue(any("not published" in line for line in logs.output))

    def test_timeout_configuration(self):
        """The datastore passes its configured timeout to the dispatcher,
        and the dispatcher defaults to `DEFAULT_TIMEOUT` when not given one.
        """
        self.assertEqual(self.butler._datastore.timeout, 5.5)
        self.assertEqual(self.butler._datastore._dispatcher.timeout, 5.5)
        self.assertEqual(SasquatchDispatcher("http://test.local", "na").timeout, DEFAULT_TIMEOUT)

    def test_requests_use_timeout(self):
        """Every HTTP call the dispatcher makes carries its timeout, so an
        unresponsive proxy cannot hang the caller indefinitely.
        """
        dispatcher = SasquatchDispatcher("http://test.local", "na", timeout=12.5)
        session = _fakeSession(getResponse=_fakeResponse(json={"data": [{"cluster_id": "abc"}]}))

        with patch(HTTP_CLIENT) as mockHttpClient:
            mockHttpClient.return_value.__enter__.return_value = session
            self.assertEqual(dispatcher.clusterId, "abc")
            self.assertTrue(dispatcher._create_topic("some.metric"))

        self.assertEqual(session.get.call_args.kwargs["timeout"], 12.5)
        self.assertEqual(session.post.call_args.kwargs["timeout"], 12.5)

    def test_cluster_id_failures_raise_dispatch_failure(self):
        """Every way the cluster id lookup can go wrong surfaces as a
        `SasquatchDispatchFailure`, which the datastore knows how to handle,
        rather than as a bare `KeyError`/`IndexError`/`JSONDecodeError`.
        """
        dispatcher = SasquatchDispatcher("http://test.local", "na")
        jsonError = requests.exceptions.JSONDecodeError("Expecting value", "<html>", 0)
        badResponses = {
            "5xx with an error document": _fakeResponse(503, json={"error_code": 50301}),
            "2xx without the data key": _fakeResponse(200, json={"error_code": 50301}),
            "2xx with no clusters": _fakeResponse(200, json={"data": []}),
            "2xx with a non-JSON body": _fakeResponse(200, jsonError=jsonError),
        }
        for description, response in badResponses.items():
            with self.subTest(description):
                with patch(HTTP_CLIENT) as mockHttpClient:
                    mockHttpClient.return_value.__enter__.return_value = _fakeSession(getResponse=response)
                    with self.assertRaises(SasquatchDispatchFailure):
                        _ = dispatcher.clusterId

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


if __name__ == "__main__":
    unittest.main()
