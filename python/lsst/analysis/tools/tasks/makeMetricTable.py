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
from astropy.table import Table
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class MakeMetricTableConnections(
    AnalysisBaseConnections,
    dimensions=("skymap",),
    defaultTemplates={"metricBundleName": "objectTableCore_metrics"},
):
    data = ct.Input(
        doc="Metric bundle to read from the butler",
        name="{metricBundleName}",
        storageClass="MetricMeasurementBundle",
        deferLoad=True,
        dimensions=(
            "skymap",
            "tract",
        ),
        multiple=True,
    )

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    metricTable = ct.Output(
        doc="A summary table of all metrics by tract.",
        name="{metricBundleName}Table",
        storageClass="ArrowAstropy",
        dimensions=("skymap",),
    )


class MakeMetricTableConfig(AnalysisBaseConfig, pipelineConnections=MakeMetricTableConnections):
    pass


class MakeMetricTableTask(AnalysisPipelineTask):
    """Turn metric bundles which are per tract into a
    summary metric table.

    TO DO: DM-44485 make sure this works for visit
    level data as well.
    """

    ConfigClass = MakeMetricTableConfig
    _DefaultName = "makeMetricTable"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Take a metric bundle for all the tracts and
        seperate it into different metrics then make it
        into a table.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
        inputRefs : `lsst.pipe.base.InputQuantizedConnection`
        outputRefs : `lsst.pipe.base.OutputQuantizedConnection`
        """

        inputs = butlerQC.get(inputRefs)
        skymap = inputs["skymap"]

        # Make a list of tracts
        tracts = []
        for data in inputRefs.data:
            tracts.append(data.dataId["tract"])

        metricBundles = []
        for inputHandle in inputs["data"]:
            metricBundles.append(inputHandle.get())

        outputs = self.run(tracts, metricBundles, skymap)
        butlerQC.put(outputs, outputRefs)

    def run(self, tracts, metricBundles, skymap):
        """Take the metric bundles and expand them out,
        add the tract corner information and then make
        a table of the information.

        Parameters
        ----------
        tracts : `list`
            A list of the tracts that the metricBundles cover
        metricBundles : `list` of
            `lsst.analysis.tools.interfaces._metricMeasurementBundle.MetricMeasurementBundle`
        skymap : `lsst.skymap`

        Returns
        -------
        tMetrics : `pipe.base.Struct` containing `astropy.table.Table`
        """

        # Make an initial dict of the columns needed
        metricsDict = {}
        tract = tracts[0]
        metricsDict["tract"] = [tract]
        corners = self.getTractCorners(skymap, tract)
        metricsDict["corners"] = [corners]

        # Add the metrics for the first tract to the dict
        for name, metrics in metricBundles[0].items():
            for metric in metrics:
                fullName = f"{name}_{metric.metric_name}"
                metricValue = metric.quantity
                metricsDict[fullName] = [metricValue]

        # Go through the rest of the tracts and check
        # if any additional columns are needed and add
        # the values to the dict

        for i, metricBundle in enumerate(metricBundles[1:]):
            tract = tracts[i + 1]
            metricsDict["tract"].append(tract)
            corners = self.getTractCorners(skymap, tract)
            metricsDict["corners"].append(corners)
            metricNames = ["tract", "corners"]
            for name, metrics in metricBundle.items():
                for metric in metrics:
                    fullName = f"{name}_{metric.metric_name}"
                    metricValue = metric.quantity
                    # Check if the metric already exists in the output
                    if fullName in metricsDict.keys():
                        metricsDict[fullName].append(metricValue)
                    else:
                        values = [np.nan] * (len(metricsDict["tract"]) - 1)
                        values.append(metricValue)
                        metricsDict[fullName] = values
                    metricNames.append(fullName)

            # If a metric that existed for a previous tract does
            # not exist for this one then add a nan
            for metricName in metricsDict.keys():
                if metricName not in metricNames:
                    metricsDict[metricName].append(np.nan)

        tMetrics = pipeBase.Struct(metricTable=Table(metricsDict))
        return tMetrics

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
