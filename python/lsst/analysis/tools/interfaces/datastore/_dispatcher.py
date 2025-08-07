# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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

__all__ = ("SasquatchDispatchPartialFailure", "SasquatchDispatchFailure", "SasquatchDispatcher")

"""Sasquatch datastore"""
import calendar
import datetime
import json
import logging
import math
import re
from collections.abc import Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, cast
from uuid import UUID, uuid4

import numpy as np
import requests
from lsst.daf.butler import DatasetRef
from lsst.resources import ResourcePath
from lsst.utils.packages import getEnvironmentPackages

from ...utils import http_client

if TYPE_CHECKING:
    from .. import MetricMeasurementBundle


log = logging.getLogger(__name__)

# Constants assocated with SasquatchDispatcher
PARTITIONS = 1
REPLICATION_FACTOR = 3

IDENTIFIER_KEYS = [
    "detector",
    "patch",
    "skymap",
    "visit",
    "tract",
    "physical_filter",
    "instrument",
    "band",
    "exposure",
    "group",
    "day_obs",
]


class SasquatchDispatchPartialFailure(RuntimeError):
    """This indicates that a Sasquatch dispatch was partially successful."""

    pass


class SasquatchDispatchFailure(RuntimeError):
    """This indicates that dispatching a
    `~lsst.analysis.tool.interface.MetricMeasurementBundle` failed.
    """

    pass


def _tag2VersionTime(productStr: str) -> tuple[str, float]:
    """Determine versions and dates from the string returned from
    getEnvironmentPackages.

    The `~lsst.utils.packages.genEnvironmentPackages` function returns the
    setup version associated with a product, along with a list of tags that
    have been added to it.

    This method splits up that return string, and determines the earliest date
    associated with the setup package version.

    Parameters
    ----------
    productStr : `str`
        The product string returned from a lookup on the result of a call to
        `~lsst.utils.packages.getEnvironmentPackages`.

    Returns
    -------
    result : `tuple` of `str`, `datetime.datetime`
        The first `str` is the version of the package, and the second is the
        datetime object associated with that released version.

    Raises
    ------
    ValueError
        Raised if there are no tags which correspond to dates.
    """
    times: list[datetime.datetime] = []
    version = productStr.split()[0]
    tags: str = re.findall("[(](.*)[)]", productStr)[0]
    for tag in tags.split():
        numDots = tag.count(".")
        numUnder = tag.count("_")
        separator = "_"
        if numDots > numUnder:
            separator = "."
        match tag.split(separator):
            # Daily tag branch.
            case ("d", year, month, day):
                dt = datetime.datetime(year=int(year), month=int(month), day=int(day))
            # Weekly tag branch.
            case ("w", year, week):
                iyear = int(year)
                iweek = int(week)
                # Use 4 as the day because releases are available starting
                # on Thursday
                dayOfWeek = 4

                # Find the first week to contain a thursday in it
                cal = calendar.Calendar()
                cal.setfirstweekday(6)
                i = 0
                for i, iterWeek in enumerate(cal.monthdatescalendar(iyear, 1)):
                    if iterWeek[dayOfWeek].month == 1:
                        break
                # Handle fromisocalendar not being able to handle week 53
                # in the case were the date was going to subtract 7 days anyway
                if i and iweek == 53:
                    i = 0
                    iweek = 52
                delta = datetime.timedelta(days=7 * i)

                # Correct for a weekly being issued in the last week of the
                # previous year, as Thursdays don't always line up evenly in
                # a week / year split.
                dt = datetime.datetime.fromisocalendar(iyear, iweek, dayOfWeek) - delta
            # Skip tags that can't be understood.
            case _:
                continue
        times.append(dt)
    if len(times) == 0:
        raise ValueError("Could not find any tags corresponding to dates")
    minTime = min(times)
    minTime.replace(tzinfo=datetime.timezone.utc)
    return version, minTime.timestamp()


@dataclass
class SasquatchDispatcher:
    """This class mediates the transfer of MetricMeasurementBundles to a
    Sasquatch http kafka proxy server.
    """

    url: str
    """Url of the Sasquatch proxy server"""

    token: str
    """Authentication token used in communicating with the proxy server"""

    namespace: str = "lsst.dm"
    """The namespace in Sasquatch in which to write the uploaded metrics"""

    def __post_init__(self) -> None:
        match ResourcePath(self.url).scheme:
            case "http" | "https":
                pass
            case _:
                raise ValueError("Proxy server must be locatable with either http or https")

        self._cluster_id: str | None = None

    @property
    def clusterId(self) -> str:
        """ClusterId of the Kafka proxy

        Notes
        -----
        The cluster Id will be fetched with a network call if it is not
        already cached.
        """
        if self._cluster_id is None:
            self._populateClusterId()
        return cast(str, self._cluster_id)

    def _populateClusterId(self) -> None:
        """Get Sasquatch kafka cluster ID."""

        headers = {"content-type": "application/json"}

        try:
            with http_client() as session:
                r = session.get(f"{self.url}/v3/clusters", headers=headers)
                cluster_id = r.json()["data"][0]["cluster_id"]
                self._cluster_id = str(cluster_id)
            r.raise_for_status()
        except requests.RequestException:
            log.error("Could not retrieve the cluster id for the specified url")
            raise SasquatchDispatchFailure("Could not retrieve the cluster id for the specified url")

    def _create_topic(self, topic_name: str) -> bool:
        """Create a kafka topic in Sasquatch.

        Parameters
        ----------
        topic_name : `str`
            The name of the kafka topic to create

        Returns
        -------
        status : `bool`
            If this does not encounter an error it will return a True success
            code, else it will return a False code.

        """

        headers = {"content-type": "application/json"}

        topic_config = {
            "topic_name": f"{self.namespace}.{topic_name}",
            "partitions_count": PARTITIONS,
            "replication_factor": REPLICATION_FACTOR,
        }

        try:
            with http_client() as session:
                r = session.post(
                    f"{self.url}/v3/clusters/{self.clusterId}/topics", json=topic_config, headers=headers
                )
            r.raise_for_status()
            log.debug("Created topic %s.%s", self.namespace, topic_name)
            return True
        except requests.HTTPError as e:
            if e.response.status_code == requests.codes.bad_request:
                log.debug("Topic %s.%s already exists.", self.namespace, topic_name)
                return True
            else:
                log.error(
                    "Unknown error occurred creating kafka topic %s %s",
                    e.response.status_code,
                    e.response.reason,
                )
                return False
        except requests.RequestException:
            log.exception("Unknown error occurred creating kafka topic")
            return False

    def _generateAvroSchema(self, metric: str, record: MutableMapping[str, Any]) -> tuple[str, bool]:
        """Infer the Avro schema from the record payload.

        Parameters
        ----------
        metric : `str`
            The name of the metric
        record : `MutableMapping`
            The prepared record for which a schema is to be generated

        Returns
        -------
        resultSchema : `str`
            A json encoded string of the resulting avro schema
        errorCode : bool
            A boolean indicating if any record fields had to be trimmed because
            a suitable schema could not be generated. True if records were
            removed, False otherwise.
        """
        schema: dict[str, Any] = {"type": "record", "namespace": self.namespace, "name": metric}

        # Record if any records needed to be trimmed
        resultsTrimmed = False

        fields = list()
        # If avro schemas cant be generated for values, they should be removed
        # from the records.
        keysToRemove: list[str] = []
        for key in record:
            value = record[key]
            avroType: Mapping[str, Any]
            if "timestamp" in key:
                avroType = {"type": "double"}
            else:
                avroType = self._python2Avro(value)
                if len(avroType) == 0:
                    continue
                if avroType.get("error_in_conversion"):
                    keysToRemove.append(key)
                    resultsTrimmed = True
                    continue
            fields.append({"name": key, **avroType})

        # remove any key that failed to have schema generated
        for key in keysToRemove:
            record.pop(key)

        schema["fields"] = fields

        return json.dumps(schema), resultsTrimmed

    def _python2Avro(self, value: Any) -> Mapping:
        """Map python type to avro schema

        Parameters
        ----------
        value : `Any`
            Any python parameter.

        Returns
        -------
        result : `Mapping`
            Return a mapping that represents an entry in an avro schema.
        """
        match value:
            case float() | np.float32() | None:
                return {"type": "float", "default": 0.0}
            case str():
                return {"type": "string", "default": ""}
            case int():
                return {"type": "int", "default": 0}
            case Sequence():
                tmp = {self._python2Avro(item)["type"] for item in value}
                if len(tmp) == 0:
                    return {}
                if len(tmp) > 1:
                    log.error(
                        "Sequence contains mixed types: %s, must be homogeneous for avro conversion "
                        "skipping record",
                        tmp,
                    )
                    return {"error_in_conversion": True}
                return {"type": "array", "items": tmp.pop()}
            case _:
                log.error("Unsupported type %s, skipping record", type(value))
                return {}

    def _handleReferencePackage(self, meta: MutableMapping, bundle: MetricMeasurementBundle) -> None:
        """Check to see if there is a reference package.

        if there is a reference package, determine the datetime associated with
        this reference package. Save the package, the version, and the date to
        the common metric fields.

        Parameters
        ----------
        meta : `MutableMapping`
            A mapping which corresponds to fields which should be encoded in
            all records.
        bundle : `MetricMeasurementBundle`
            The bundled metrics
        """
        package_version, package_timestamp = "", 0.0
        if ref_package := getattr(bundle, "reference_package", ""):
            ref_package = bundle.reference_package
            packages = getEnvironmentPackages(True)
            if package_info := packages.get(ref_package):
                try:
                    package_version, package_timestamp = _tag2VersionTime(package_info)
                except ValueError:
                    # Could not extract package timestamp leaving empty
                    pass
        # explicit handle if None was set in the bundle for the package
        meta["reference_package"] = ref_package or ""
        meta["reference_package_version"] = package_version
        meta["reference_package_timestamp"] = package_timestamp

    def _handleTimes(self, meta: MutableMapping, bundle: MetricMeasurementBundle, run: str) -> None:
        """Add times to the meta fields mapping.

        Add all appropriate timestamp fields to the meta field mapping. These
        will be added to all records.

        This method will also look at the bundle to see if it defines a
        preferred time. It so it sets that time as the main time stamp to be
        used for this record.

        Parameters
        ----------
        meta : `MutableMapping`
            A mapping which corresponds to fields which should be encoded in
            all records.
        bundle : `MetricMeasurementBundle`
            The bundled metrics
        run : `str`
            The `~lsst.daf.butler.Butler` collection where the
            `MetricMeasurementBundle` is stored.
        """
        # Determine the timestamp associated with the run, if someone abused
        # the run collection, use the current timestamp
        if re.match(r"\d{8}T\d{6}Z", stamp := run.split("/")[-1]):
            run_timestamp = datetime.datetime.strptime(stamp, r"%Y%m%dT%H%M%S%z")
        else:
            run_timestamp = datetime.datetime.now()
        meta["run_timestamp"] = run_timestamp.timestamp()

        # If the bundle supports supplying timestamps, dispatch on the type
        # specified.
        if hasattr(bundle, "timestamp_version") and bundle.timestamp_version:
            match bundle.timestamp_version:
                case "reference_package_timestamp":
                    if not meta["reference_package_timestamp"]:
                        log.error("Reference package timestamp is empty, using run_timestamp")
                        meta["timestamp"] = meta["run_timestamp"]
                    else:
                        meta["timestamp"] = meta["reference_package_timestamp"]
                case "run_timestamp":
                    meta["timestamp"] = meta["run_timestamp"]
                case "current_timestamp":
                    timeStamp = datetime.datetime.now()
                    meta["timestamp"] = timeStamp.timestamp()
                case "dataset_timestamp":
                    log.error("dataset timestamps are not yet supported, run_timestamp will be used")
                    meta["timestamp"] = meta["run_timestamp"]
                case str(value) if "explicit_timestamp" in value:
                    try:
                        _, splitTime = value.split(":")
                    except ValueError as excpt:
                        raise ValueError(
                            "Explicit timestamp must be given in the format 'explicit_timestamp:datetime', "
                            "where datetime is given in the form '%Y%m%dT%H%M%S%z"
                        ) from excpt
                    meta["timestamp"] = datetime.datetime.strptime(splitTime, r"%Y%m%dT%H%M%S%z").timestamp()
                case _:
                    log.error(
                        "Timestamp version %s is not supported, run_timestamp will be used",
                        bundle.timestamp_version,
                    )
                    meta["timestamp"] = meta["run_timestamp"]
        # Default to using the run_timestamp.
        else:
            meta["timestamp"] = meta["run_timestamp"]

    def _handleIdentifier(
        self,
        meta: MutableMapping,
        identifierFields: Mapping[str, Any] | None,
        datasetIdentifier: str | None,
        bundle: MetricMeasurementBundle,
    ) -> None:
        """Add an identifier to the meta record mapping.

        If the bundle declares a dataset identifier to use add that to the
        record, otherwise use 'Generic' as the identifier. If the
        datasetIdentifier parameter is specified, that is used instead of
        anything specified by the bundle.

        This will also add any identifier fields supplied to the meta record
        mapping.

        Together these values (in addition to the timestamp and topic) should
        uniquely identify an upload to the Sasquatch system.

        Parameters
        ----------
        meta : `MutableMapping`
            A mapping which corresponds to fields which should be encoded in
            all records.
        identifierFields: `Mapping` or `None`
            The keys and values in this mapping will be both added as fields
            in the record, and used in creating a unique tag for the uploaded
            dataset type. I.e. the timestamp, and the tag will be unique, and
            each record will belong to one combination of such.
        datasetIdentifier : `str` or `None`
            A string which will be used in creating unique identifier tags.
        bundle : `MetricMeasurementBundle`
            The bundle containing metric values to upload.
        """
        identifier: str
        if datasetIdentifier is not None:
            identifier = datasetIdentifier
        elif hasattr(bundle, "dataset_identifier") and bundle.dataset_identifier is not None:
            identifier = bundle.dataset_identifier
        else:
            identifier = "Generic"

        meta["dataset_tag"] = identifier

        if identifierFields is None:
            identifierFields = {}
        for key in IDENTIFIER_KEYS:
            value = identifierFields.get(key, "")
            meta[key] = f"{value}"

    def _prepareBundle(
        self,
        bundle: MetricMeasurementBundle,
        run: str,
        datasetType: str,
        timestamp: datetime.datetime | None = None,
        id: UUID | None = None,
        identifierFields: Mapping | None = None,
        datasetIdentifier: str | None = None,
        extraFields: Mapping | None = None,
    ) -> tuple[Mapping[str, list[Any]], bool]:
        """Encode all of the inputs into a format that can be sent to the
        kafka proxy server.

        Parameters
        ----------
        bundle : `MetricMeasurementBundle`
            The bundle containing metric values to upload.
        run : `str`
            The run name to associate with these metric values. If this bundle
            is also stored in the butler, this should be the butler run
            collection the bundle is stored in the butler.
        datasetType : `str`
            The dataset type name associated with this
            `MetricMeasurementBundle`
        timestamp : `datetime.datetime`, optional
            The timestamp to be associated with the measurements in the ingress
            database. If this value is None, timestamp will be set by the run
            time or current time.
        id : `UUID`, optional
            The UUID of the `MetricMeasurementBundle` within the butler. If
            `None`, a new random UUID will be generated so that each record in
            Sasquatch will have a unique value.
        identifierFields: `Mapping`, optional
            The keys and values in this mapping will be both added as fields
            in the record, and used in creating a unique tag for the uploaded
            dataset type. I.e. the timestamp, and the tag will be unique, and
            each record will belong to one combination of such.
        datasetIdentifier : `str`, optional
            A string which will be used in creating unique identifier tags.
        extraFields: `Mapping`, optional
            Extra mapping keys and values that will be added as fields to the
            dispatched record.

        Returns
        -------
        result : `Mapping` of `str` to `list`
            A mapping of metric name of list of metric measurement records.
        status : `bool`
            A status boolean indicating if some records had to be skipped due
            to a problem parsing the bundle.
        """
        if id is None:
            id = uuid4()
        sid = str(id)
        meta: dict[str, Any] = dict()

        # Add other associated common fields
        meta["id"] = sid
        meta["run"] = run
        meta["dataset_type"] = datasetType

        # Check to see if the bundle declares a reference package
        self._handleReferencePackage(meta, bundle)

        # Handle the various timestamps that could be associated with a record
        self._handleTimes(meta, bundle, run)

        # Always use the supplied timestamp if one was passed to use.
        if timestamp is not None:
            meta["timestamp"] = timestamp.timestamp()

        self._handleIdentifier(meta, identifierFields, datasetIdentifier, bundle)

        # Add in any other fields that were supplied to the function call.
        if extraFields is not None:
            meta.update(extraFields)

        metricRecords: dict[str, list[Any]] = dict()

        # Record if any records needed skipped
        resultsTrimmed = False

        # Look at each of the metrics in the bundle (name, values)
        for metric, measurements in bundle.items():
            # Create a list which will contain the records for each measurement
            # associated with metric.
            metricRecordList = metricRecords.setdefault(f"{bundle.metricNamePrefix}{metric}", list())

            record: dict[str, Any] = meta.copy()

            # loop over each metric measurement within the metric
            for measurement in measurements:
                # need to extract any tags, package info, etc
                note_key = f"{measurement.metric_name.metric}.metric_tags"
                record["tags"] = dict(measurement.notes.items()).get(note_key, list())

                # Missing values are replaced by 0 in sasquatch, see RFC-763.
                name = ""
                value = 0.0
                match measurement.json:
                    case {"metric": name, "value": None}:
                        pass
                    case {"metric": name, "value": value}:
                        if math.isnan(value):
                            log.error(
                                "Measurement %s had a value that is a NaN, dispatch will be skipped",
                                measurement,
                            )
                            resultsTrimmed = True
                            continue
                        # JSON will not serialize np.float32, must cast.
                        if isinstance(value, np.float32):
                            value = float(value)
                        pass
                    case {"value": _}:
                        log.error("Measurement %s does not contain the key 'metric'", measurement)
                        resultsTrimmed = True
                        continue
                    case {"metric": _}:
                        log.error("Measurement %s does not contain the key 'value'", measurement)
                        resultsTrimmed = True
                        continue
                record[name] = value

            metricRecordList.append({"value": record})
        return metricRecords, resultsTrimmed

    def dispatch(
        self,
        bundle: MetricMeasurementBundle,
        run: str,
        datasetType: str,
        timestamp: datetime.datetime | None = None,
        id: UUID | None = None,
        datasetIdentifier: str | None = None,
        identifierFields: Mapping | None = None,
        extraFields: Mapping | None = None,
    ) -> None:
        """Dispatch a `MetricMeasurementBundle` to Sasquatch.

        Parameters
        ----------
        bundle : `MetricMeasurementBundle`
            The bundle containing metric values to upload.
        run : `str`
            The run name to associate with these metric values. If this bundle
            is also stored in the butler, this should be the butler run
            collection the bundle is stored in the butler. This will be used
            in generating uniqueness constraints in Sasquatch.
        datasetType : `str`
            The dataset type name associated with this
            `MetricMeasurementBundle`.
        timestamp : `datetime.datetime`, optional
            The timestamp to be associated with the measurements in the ingress
            database. If this value is None, timestamp will be set by the run
            time or current time.
        id : `UUID`, optional
            The UUID of the `MetricMeasurementBundle` within the Butler. If
            `None`, a new random UUID will be generated so that each record in
            Sasquatch will have a unique value.
        datasetIdentifier : `str`, optional
            A string which will be used in creating unique identifier tags. If
            `None`, a default value will be inserted.
        identifierFields: `Mapping`, optional
            The keys and values in this mapping will be both added as fields
            in the record, and used in creating a unique tag for the uploaded
            dataset type. I.e. the timestamp, and the tag will be unique, and
            each record will belong to one combination of such. Examples of
            entries would be things like visit or tract.
        extraFields: `Mapping`, optional
            Extra mapping keys and values that will be added as fields to the
            dispatched record.

        Raises
        ------
        SasquatchDispatchPartialFailure, SasquatchDispatchFailure
            Raised if there were any errors in dispatching a bundle.
        """
        if id is None:
            id = uuid4()

        # Prepare the bundle by transforming it to a list of metric records
        metricRecords, recordsTrimmed = self._prepareBundle(
            bundle=bundle,
            run=run,
            datasetType=datasetType,
            timestamp=timestamp,
            id=id,
            datasetIdentifier=datasetIdentifier,
            identifierFields=identifierFields,
            extraFields=extraFields,
        )

        headers = {"content-type": "application/vnd.kafka.avro.v2+json"}
        data: dict[str, Any] = dict()
        partialUpload = False
        uploadFailed = []

        with http_client() as session:
            for metric, record in metricRecords.items():
                # create the kafka topic if it does not already exist
                if not self._create_topic(metric):
                    log.error("Topic not created, skipping dispatch of %s", metric)
                    continue
                recordValue = record[0]["value"]
                # Generate schemas for each record
                data["value_schema"], schemaTrimmed = self._generateAvroSchema(metric, recordValue)
                data["records"] = record

                if schemaTrimmed:
                    partialUpload = True

                try:
                    r = session.post(
                        f"{self.url}/topics/{self.namespace}.{metric}", json=data, headers=headers
                    )
                    r.raise_for_status()
                    log.debug("Succesfully sent data for metric %s", metric)
                    uploadFailed.append(False)
                except requests.HTTPError as e:
                    log.error(
                        "There was a problem submitting the metric %s: %s, %s",
                        metric,
                        e.response.status_code,
                        e.response.reason,
                    )
                    uploadFailed.append(True)
                    partialUpload = True
                except requests.RequestException as e:
                    # Don't log full stack trace because there may be lots
                    # of these.
                    log.error("There was a problem submitting the metric %s: %s", metric, e)
                    uploadFailed.append(True)
                    partialUpload = True

        # There may be no metrics to try to upload, and thus the uploadFailed
        # list may be empty, check before issuing failure
        if len(uploadFailed) > 0 and all(uploadFailed):
            raise SasquatchDispatchFailure("All records were unable to be uploaded.")

        if partialUpload or recordsTrimmed:
            raise SasquatchDispatchPartialFailure("One or more records may not have been uploaded entirely")

    def dispatchRef(
        self,
        bundle: MetricMeasurementBundle,
        ref: DatasetRef,
        timestamp: datetime.datetime | None = None,
        extraFields: Mapping | None = None,
        datasetIdentifier: str | None = None,
    ) -> None:
        """Dispatch a `MetricMeasurementBundle` to Sasquatch with a known
        `DatasetRef`.

        Parameters
        ----------
        bundle : `MetricMeasurementBundle`
            The bundle containing metric values to upload.
        ref : `DatasetRef`
            The `Butler` dataset ref corresponding to the input
            `MetricMeasurementBundle`.
        timestamp : `datetime.datetime`, optional
            The timestamp to be associated with the measurements in the ingress
            database. If this value is None, timestamp will be set by the run
            time or current time.
        extraFields: `Mapping`, optional
            Extra mapping keys and values that will be added as fields to the
            dispatched record if not None.
        datasetIdentifier : `str`, optional
            A string which will be used in creating unique identifier tags. If
            None, a default value will be inserted.

        Raises
        ------
        SasquatchDispatchPartialFailure, SasquatchDispatchFailure
            Raised if there were any errors in dispatching a bundle.
        """
        # Parse the relevant info out of the dataset ref.
        serializedRef = ref.to_simple()
        id = serializedRef.id
        if serializedRef.run is None:
            run = "<unknown>"
        else:
            run = serializedRef.run
        dstype = serializedRef.datasetType
        datasetType = dstype.name if dstype is not None else ""
        dataRefMapping = serializedRef.dataId.dataId if serializedRef.dataId else None

        self.dispatch(
            bundle,
            run=run,
            timestamp=timestamp,
            datasetType=datasetType,
            id=id,
            identifierFields=dataRefMapping,
            extraFields=extraFields,
            datasetIdentifier=datasetIdentifier,
        )
