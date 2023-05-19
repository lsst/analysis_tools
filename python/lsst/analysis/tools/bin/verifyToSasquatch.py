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

__all__ = [
    "main",
]

import argparse
from collections import defaultdict
from collections.abc import Iterable, Mapping

import lsst.verify
from lsst.analysis.tools.interfaces import MetricMeasurementBundle
from lsst.daf.butler import Butler, DataCoordinate, DatasetRef


def makeParser():
    parser = argparse.ArgumentParser()
    return parser


def _bundle_metrics(
    butler: Butler, metricValues: Iterable[DatasetRef]
) -> Mapping[tuple[str, str, DataCoordinate], MetricMeasurementBundle]:
    """Organize free metric values into metric bundles while preserving as much
    information as practical.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The Butler repository containing the metric values.
    metricValues : `~collections.abc.Iterable` [`lsst.daf.butler.DatasetRef`]
        The datasets to bundle. All references must point to ``MetricValue``
        datasets.

    Returns
    -------
    bundles : `~collections.abc.Mapping`
        A collection of
        `lsst.analysis.tools.interfaces.MetricMeasurementBundle`, one for each
        combination of distinct metadata. The mapping key is a tuple of (run,
        dataset type, data ID), and the value is the corresponding bundle.
        To simplify the uploaded schemas, the bundle uses metrics' relative
        (unqualified) names even if the input measurements were
        fully-qualified.
    """
    bundles = defaultdict(MetricMeasurementBundle)
    for ref in metricValues:
        value = butler.get(ref)
        # MeasurementMetricBundle doesn't validate input.
        if not isinstance(value, lsst.verify.Measurement):
            raise ValueError(f"{ref} is not a metric value.")

        # HACK: in general, metric names are fully qualified, and this becomes
        # the InfluxDB field name. lsst.verify-style metrics have unique names
        # already, so remove the package qualification.
        value = lsst.verify.Measurement(
            value.metric_name.metric, value.quantity, value.blobs.values(), value.extras, value.notes
        )
        # These metrics weren't created by actions. Sasquatch requires that
        # each actionId produce the same metrics on every run (see
        # https://sasquatch.lsst.io/user-guide/avro.html), so choose something
        # unique to the metric.
        actionId = value.metric_name.metric

        bundle = bundles[(ref.run, ref.datasetType.name, ref.dataId)]
        bundle.setdefault(actionId, []).append(value)
    return bundles


def main():
    makeParser().parse_args()
