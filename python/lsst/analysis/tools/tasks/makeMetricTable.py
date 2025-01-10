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
from __future__ import annotations

__all__ = (
    "MakeMetricTableConfig",
    "MakeMetricTableTask",
)

import lsst.pipe.base as pipeBase
import numpy as np
from astropy import units as u
from astropy.table import Table
from lsst.pex.config import ListField
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap


class MakeMetricTableConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=(),
    defaultTemplates={"metricBundleName": "", "outputTableName": ""},
):
    data = ct.Input(
        doc="Metric bundle to read from the butler",
        name="{metricBundleName}",
        storageClass="MetricMeasurementBundle",
        deferLoad=True,
        dimensions=(),
        multiple=True,
    )

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    metricTable = ct.Output(
        doc="Table containing metrics, one row per input metric bundle",
        name="{outputTableName}",
        storageClass="ArrowAstropy",
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        self.data = ct.Input(
            doc=self.data.doc,
            name=self.data.name,
            storageClass=self.data.storageClass,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputDataDimensions)),
            multiple=self.data.multiple,
        )
        self.metricTable = ct.Output(
            doc=self.metricTable.doc,
            name=self.metricTable.name,
            storageClass=self.metricTable.storageClass,
            dimensions=frozenset(sorted(config.outputTableDimensions)),
        )

        assert config is not None, "Missing required config object."

        if "tract" not in config.inputDataDimensions:
            del self.skymap

        self.dimensions.update(frozenset(sorted(config.outputTableDimensions)))


class MakeMetricTableConfig(
    pipeBase.PipelineTaskConfig,
    pipelineConnections=MakeMetricTableConnections,
):
    inputDataDimensions = ListField[str](
        doc="Dimensions of the input data.",
        default=(),
        optional=False,
    )
    outputTableDimensions = ListField[str](
        doc="Dimensions of the output data.",
        default=(),
        optional=False,
    )
    dataIdFieldsToIncludeAsColumns = ListField[str](
        doc="DataId fields to include as columns in the table. "
        "These are added in addition to the Metric names. "
        "At least one field must be specified.",
        default=None,
        optional=False,
    )


class MakeMetricTableTask(pipeBase.PipelineTask):
    """Turn metric bundles and combine them into a metric table."""

    ConfigClass = MakeMetricTableConfig
    _DefaultName = "makeMetricTable"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Take a set of metric bundles, seperate each into its different
        metrics, then put the values into a table with the metric names as
        column headers.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
        inputRefs : `lsst.pipe.base.InputQuantizedConnection`
        outputRefs : `lsst.pipe.base.OutputQuantizedConnection`
        """

        inputs = butlerQC.get(inputRefs)
        if "skymap" in inputs:
            skymap = inputs["skymap"]
        else:
            skymap = None

        # Extract the info from the dataIds that is needed
        # to populate the requested columns.
        fields = self.config.dataIdFieldsToIncludeAsColumns
        dataIdInfo = []
        for data in inputRefs.data:
            dataIdInfo.append({field: data.dataId[field] for field in fields})

        metricBundles = []
        for inputHandle in inputs["data"]:
            metricBundles.append(inputHandle.get())

        outputs = self.run(dataIdInfo, metricBundles, skymap)
        butlerQC.put(outputs, outputRefs)

    def run(self, dataIdInfo, metricBundles, skymap):
        """Take the metric bundles and expand them out, then make a table of
        the information. Add tract corner information if the bundles are
        tract-level.

        Parameters
        ----------
        dataIdInfo : `list`
            A list of dicts that hold information extracted from the metric
            bundle dataIds.
        metricBundles : `list` of
            `lsst.analysis.tools.interfaces._metricMeasurementBundle.MetricMeasurementBundle`
        skymap : `lsst.skymap`

        Returns
        -------
        metricTableStruct : `pipe.base.Struct` containing `astropy.table.Table`
        """

        if len(dataIdInfo) == 0:
            raise pipeBase.NoWorkFound("dataIdInfo list is empty")
        if len(metricBundles) == 0:
            raise pipeBase.NoWorkFound("metricBundles list is empty")

        metricsDict = {}
        metricUnits = {}

        # Add requested info from the first dataId to the metrics dict.
        for key, value in dataIdInfo[0].items():
            metricsDict[key] = [value]

        # Add tract corners if inputs are at the tract-level
        if "tract" in self.config.inputDataDimensions:
            corners = self.getTractCorners(skymap, dataIdInfo[0]["tract"])
            metricsDict["corners"] = [corners]

        # Add the metrics and units from the first bundle to the dicts
        for name, metrics in metricBundles[0].items():
            percentUnitUsed = False
            for metric in metrics:
                fullName = f"{name}_{metric.metric_name}"
                metricsDict[fullName] = [metric.quantity.value]
                # "Dimensionless" and percent units not allowed in Tables:
                if metric.quantity.unit is u.dimensionless_unscaled:
                    continue
                elif metric.quantity.unit is u.pct:
                    percentUnitUsed = True
                    self.log.debug(
                        "Unable to propagate astropy percent unit for metric %s "
                        "in bundle %s to astropy table.",
                        metric.metric_name,
                        name,
                    )
                    continue
                else:
                    metricUnits[fullName] = metric.quantity.unit
            if percentUnitUsed:
                self.log.warn(
                    "One or more metrics in the %s metric bundle uses the percent unit, "
                    "which is not supported in astropy Tables. "
                    "The value(s) have been propagated without the percent unit. "
                    "Use --log-level debug to list all affected metrics.",
                    name,
                )

        # Check if any additional columns are needed; add to dict if needed.
        for i, metricBundle in enumerate(metricBundles[1:]):
            for key, value in dataIdInfo[i + 1].items():
                metricsDict[key].append(value)

            if "tract" in self.config.inputDataDimensions:
                corners = self.getTractCorners(skymap, dataIdInfo[0]["tract"])
                metricsDict["corners"].append(corners)

            metricNames = list(metricsDict)
            for name, metrics in metricBundle.items():
                for metric in metrics:
                    fullName = f"{name}_{metric.metric_name}"
                    # Check if the metric already exists in the output
                    if fullName in metricsDict.keys():
                        metricsDict[fullName].append(metric.quantity.value)
                    else:
                        values = [np.nan] * (len(metricsDict[metricNames[0]]) - 1)
                        values.append(metric.quantity.value)
                        if (
                            metric.quantity.unit is not u.dimensionless_unscaled
                            and metric.quantity.unit is not u.pct
                        ):
                            metricUnits[fullName] = metric.quantity.unit
                        metricsDict[fullName] = values
                    metricNames.append(fullName)

            # If a metric that existed in a previous bundle does
            # not exist for this one then add a nan
            for metricName in metricsDict.keys():
                if metricName not in metricNames:
                    metricsDict[metricName].append(np.nan)

        metricTableStruct = pipeBase.Struct(metricTable=Table(metricsDict, units=metricUnits))
        return metricTableStruct

    def getTractCorners(self, skymap, tract):
        """Calculate the corners of a tract, given  skymap.

        Parameters
        ----------
        skymap : `lsst.skymap`
        tract : `int`

        Returns
        -------
        corners : `list` of `tuples` of `float`

        Notes
        -----
        Corners are returned in degrees and wrapped in ra.
        """
        # Find the tract corners
        tractCorners = skymap[tract].getVertexList()
        corners = [(corner.getRa().asDegrees(), corner.getDec().asDegrees()) for corner in tractCorners]
        minRa = np.min([corner[0] for corner in corners])
        maxRa = np.max([corner[0] for corner in corners])
        # If the tract needs wrapping in ra, wrap it
        if maxRa - minRa > 10:
            x = maxRa
            maxRa = 360 + minRa
            minRa = x
            minDec = np.min([corner[1] for corner in corners])
            maxDec = np.max([corner[1] for corner in corners])
            corners = [(minRa, minDec), (maxRa, minDec), (maxRa, maxDec), (minRa, maxDec)]

        return corners
