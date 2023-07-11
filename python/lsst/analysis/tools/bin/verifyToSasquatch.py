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
import copy
import datetime
import logging
from collections import defaultdict
from collections.abc import Iterable, Mapping

import lsst.verify
from lsst.analysis.tools.interfaces import MetricMeasurementBundle
from lsst.analysis.tools.interfaces.datastore import SasquatchDispatcher
from lsst.daf.butler import Butler, DataCoordinate, DatasetRef

logging.basicConfig()
_LOG = logging.getLogger(__name__)
_LOG.setLevel(logging.INFO)


_BASE_URL = "https://usdf-rsp-dev.slac.stanford.edu/sasquatch-rest-proxy/"


def makeParser():
    parser = argparse.ArgumentParser(
        description="""Upload Measurement datasets from a Butler repository to Sasquatch.

                    This script handles metric values persisted directly using
                    lsst.verify tooling. It is neither necessary nor useful for
                    MetricMeasurementBundles created using analysis_tools
                    tooling, and is provided solely for backwards compatibility
                    with the older system.
                    """,
        add_help=True,
    )
    parser.add_argument("repo", help="The Butler repository from which to upload metric values.")
    parser.add_argument(
        "collections",
        action="extend",
        nargs="+",
        help="The collection(s) in which to search for metric values. These can "
        "be specified in any notation recognized by Middleware.",
    )
    parser.add_argument("--dataset", required=True, help="The dataset on which the metrics were measured.")
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run this command while uploading to the lsst.debug test "
        "namespace. Any --namespace argument is ignored.",
    )
    parser.add_argument(
        "--where", default="", help="Butler query to select metric values for upload (default: all values)."
    )
    parser.add_argument(
        "--date-created",
        type=datetime.datetime.fromisoformat,
        help="ISO8601 formatted datetime in UTC for the Measurement creation "
        "date, e.g. 2021-06-30T22:28:25Z. If not provided, the run time or "
        "current time is used.",
    )

    api_group = parser.add_argument_group("Sasquatch API arguments")
    api_group.add_argument(
        "--namespace",
        default="lsst.dm",
        help="The Sasquatch namespace to which to upload the metric values (default: lsst.dm)",
    )
    api_group.add_argument(
        "--url",
        dest="base_url",
        default=_BASE_URL,
        help=f"Root URL of Sasquatch proxy server (default: {_BASE_URL}).",
    )
    api_group.add_argument("--token", default="na", help="Authentication token for the proxy server.")

    return parser


class _AppendDict(argparse.Action):
    """An action analogous to the build-in 'append' that appends to a `dict`
    instead of a `list`.

    Inputs are assumed to be strings in the form "key=value"; any input that
    does not contain exactly one "=" character is invalid. If the default value
    is non-empty, the default key-value pairs may be overwritten by values from
    the command line.
    """

    def __init__(
        self,
        option_strings,
        dest,
        nargs=None,
        const=None,
        default=None,
        type=None,
        choices=None,
        required=False,
        help=None,
        metavar=None,
    ):
        if default is None:
            default = {}
        if not isinstance(default, Mapping):
            argname = option_strings if option_strings else metavar if metavar else dest
            raise TypeError(f"Default for {argname} must be a mapping or None, got {default!r}.")
        super().__init__(option_strings, dest, nargs, const, default, type, choices, required, help, metavar)

    def __call__(self, parser, namespace, values, option_string=None):
        # argparse doesn't make defensive copies, so namespace.dest may be
        # the same object as self.default. Do the copy ourselves and avoid
        # modifying the object previously in namespace.dest.
        mapping = copy.copy(getattr(namespace, self.dest))

        # Sometimes values is a copy of default instead of an input???
        if isinstance(values, Mapping):
            mapping.update(values)
        else:
            # values may be either a string or list of strings, depending on
            # nargs. Unsafe to test for Sequence, because a scalar string
            # passes.
            if not isinstance(values, list):
                values = [values]
            for value in values:
                vars = value.split("=")
                if len(vars) != 2:
                    raise ValueError(f"Argument {value!r} does not match format 'key=value'.")
                mapping[vars[0]] = vars[1]

        # Other half of the defensive copy.
        setattr(namespace, self.dest, mapping)


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
    args = makeParser().parse_args()
    if args.test:
        args.namespace = "lsst.debug"

    butler = Butler(args.repo, collections=args.collections, writeable=False)
    metricTypes = {t for t in butler.registry.queryDatasetTypes() if t.storageClass_name == "MetricValue"}
    metricValues = butler.registry.queryDatasets(metricTypes, where=args.where, findFirst=True)
    _LOG.info("Found %d metric values in %s.", metricValues.count(), args.collections)

    bundles = _bundle_metrics(butler, metricValues)
    dispatcher = SasquatchDispatcher(url=args.base_url, token=args.token, namespace=args.namespace)
    _LOG.info("Uploading to %s @ %s...", args.namespace, args.base_url)
    for (run, datasetType, dataId), bundle in bundles.items():
        dispatcher.dispatch(
            bundle,
            run=run,
            datasetType=datasetType,
            timestamp=args.date_created,
            datasetIdentifier=args.dataset,
            identifierFields=dataId,
        )
