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

import tempfile
import unittest

import astropy.units as u

import lsst.daf.butler.tests as butlerTests
from lsst.analysis.tools.bin.verifyToSasquatch import _bundle_metrics
from lsst.analysis.tools.interfaces import MetricMeasurementBundle
from lsst.daf.butler import CollectionType, DataCoordinate
from lsst.verify import Measurement


class VerifyToSasquatchTestSuite(unittest.TestCase):
    def setUp(self):
        super().setUp()

        repo = tempfile.TemporaryDirectory()
        # TemporaryDirectory warns on leaks; addCleanup also keeps the TD from
        # getting garbage-collected.
        self.addCleanup(tempfile.TemporaryDirectory.cleanup, repo)
        self.butler = butlerTests.makeTestRepo(repo.name)

        butlerTests.addDataIdValue(self.butler, "instrument", "notACam")
        for visit in [42, 43]:
            butlerTests.addDataIdValue(self.butler, "visit", visit)
            for detector in range(4):
                butlerTests.addDataIdValue(self.butler, "detector", detector)
        butlerTests.addDatasetType(
            self.butler,
            "metricvalue_nopackage_fancyMetric",
            {"instrument", "visit", "detector"},
            "MetricValue",
        )
        butlerTests.addDatasetType(
            self.butler,
            "metricvalue_nopackage_fancierMetric",
            {"instrument", "visit", "detector"},
            "MetricValue",
        )
        butlerTests.addDatasetType(
            self.butler, "metricvalue_nopackage_plainMetric", {"instrument"}, "MetricValue"
        )
        self.butler.registry.registerCollection("run1", CollectionType.RUN)
        self.butler.registry.registerCollection("run2", CollectionType.RUN)

    def _standardize(self, dataId):
        """Convert an arbitrary data ID to a DataCoordinate."""
        return DataCoordinate.standardize(dataId, universe=self.butler.dimensions)

    def test_bundle_metrics_nometrics(self):
        refs = self.butler.registry.queryDatasets("metricvalue_nopackage_*", collections=...)
        bundles = _bundle_metrics(self.butler, refs)
        # MetricMeasurementBundle.equals ignores metadata
        self.assertDictEqual(bundles, {})

    def test_bundle_metrics_onemetric(self):
        m = Measurement("nopackage.fancyMetric", 42.2 * u.s)
        self.butler.put(
            m, "metricvalue_nopackage_fancyMetric", run="run1", instrument="notACam", visit=42, detector=2
        )
        refs = self.butler.registry.queryDatasets("metricvalue_nopackage_fancyMetric", collections=...)
        bundles = _bundle_metrics(self.butler, refs)
        # MetricMeasurementBundle.equals ignores metadata
        self.assertDictEqual(
            bundles,
            {
                (
                    "run1",
                    "metricvalue_nopackage_fancyMetric",
                    self._standardize({"instrument": "notACam", "visit": 42, "detector": 2}),
                ): MetricMeasurementBundle({"fancyMetric": [Measurement(m.metric_name.metric, m.quantity)]}),
            },
        )

    def test_bundle_metrics_onemetrictwice(self):
        m2 = Measurement("nopackage.fancyMetric", 42.2 * u.s)
        self.butler.put(
            m2, "metricvalue_nopackage_fancyMetric", run="run1", instrument="notACam", visit=42, detector=2
        )
        m1 = Measurement("nopackage.fancyMetric", 43.1 * u.m)
        self.butler.put(
            m1, "metricvalue_nopackage_fancyMetric", run="run1", instrument="notACam", visit=43, detector=1
        )
        refs = self.butler.registry.queryDatasets("metricvalue_nopackage_fancyMetric", collections=...)
        bundles = _bundle_metrics(self.butler, refs)
        # MetricMeasurementBundle.equals ignores metadata
        self.assertDictEqual(
            bundles,
            {
                (
                    "run1",
                    "metricvalue_nopackage_fancyMetric",
                    self._standardize({"instrument": "notACam", "visit": 42, "detector": 2}),
                ): MetricMeasurementBundle(
                    {"fancyMetric": [Measurement(m2.metric_name.metric, m2.quantity)]}
                ),
                (
                    "run1",
                    "metricvalue_nopackage_fancyMetric",
                    self._standardize({"instrument": "notACam", "visit": 43, "detector": 1}),
                ): MetricMeasurementBundle(
                    {"fancyMetric": [Measurement(m1.metric_name.metric, m1.quantity)]}
                ),
            },
        )

    def test_bundle_metrics_threemetrics(self):
        m2 = Measurement("nopackage.fancyMetric", 42.2 * u.s)
        self.butler.put(
            m2, "metricvalue_nopackage_fancyMetric", run="run1", instrument="notACam", visit=42, detector=2
        )
        m1 = Measurement("nopackage.fancierMetric", 43.1 * u.m)
        self.butler.put(
            m1, "metricvalue_nopackage_fancierMetric", run="run2", instrument="notACam", visit=43, detector=1
        )
        m_ = Measurement("nopackage.plainMetric", 127.0 * u.kg)
        self.butler.put(m_, "metricvalue_nopackage_plainMetric", run="run1", instrument="notACam")
        refs = self.butler.registry.queryDatasets("metricvalue_nopackage_*", collections=...)
        bundles = _bundle_metrics(self.butler, refs)
        # MetricMeasurementBundle.equals ignores metadata
        self.assertDictEqual(
            bundles,
            {
                (
                    "run1",
                    "metricvalue_nopackage_fancyMetric",
                    self._standardize({"instrument": "notACam", "visit": 42, "detector": 2}),
                ): MetricMeasurementBundle(
                    {"fancyMetric": [Measurement(m2.metric_name.metric, m2.quantity)]}
                ),
                (
                    "run2",
                    "metricvalue_nopackage_fancierMetric",
                    self._standardize({"instrument": "notACam", "visit": 43, "detector": 1}),
                ): MetricMeasurementBundle(
                    {"fancierMetric": [Measurement(m1.metric_name.metric, m1.quantity)]}
                ),
                (
                    "run1",
                    "metricvalue_nopackage_plainMetric",
                    self._standardize({"instrument": "notACam"}),
                ): MetricMeasurementBundle(
                    {"plainMetric": [Measurement(m_.metric_name.metric, m_.quantity)]}
                ),
            },
        )

    def test_bundle_metrics_badmetric(self):
        butlerTests.addDatasetType(
            self.butler, "metricvalue_nopackage_notAMetric", {"instrument", "visit"}, "StructuredDataDict"
        )
        self.butler.put(
            {"foo": "bar"}, "metricvalue_nopackage_notAMetric", run="run1", instrument="notACam", visit=42
        )
        refs = self.butler.registry.queryDatasets("metricvalue_nopackage_notAMetric", collections=...)
        with self.assertRaises(ValueError):
            _bundle_metrics(self.butler, refs)
