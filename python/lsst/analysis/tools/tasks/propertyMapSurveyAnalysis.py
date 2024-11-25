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

__all__ = [
    "PropertyMapSurveyWideAnalysisConfig",
    "PropertyMapSurveyWideAnalysisTask",
]

from typing import Any, Mapping, Union

from lsst.daf.butler import DataCoordinate
from lsst.pipe.base import InputQuantizedConnection, OutputQuantizedConnection, QuantumContext
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class PropertyMapSurveyWideAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "band"),
    defaultTemplates={"outputName": "propertyMapSurvey"},
):
    healSparsePropertyMapsConfig = ct.Input(
        doc="Configuration parameters for HealSparseInputMapTask in pipe_tasks.",
        name="healSparsePropertyMaps_config",
        storageClass="Config",
        dimensions=(),
    )

    skymap = ct.Input(
        doc="The skymap which the data has been mapped onto.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        operationNameLookup = {
            "min": "Minimum",
            "max": "Maximum",
            "mean": "Mean",
            "weighted_mean": "Weighted mean",
            "sum": "Sum",
        }

        # Making connections for the maps that are configured to run.
        for name in config.atools.fieldNames:
            propertyName, operationName = name.split("_consolidated_map_")
            coaddName, propertyName = propertyName.split("Coadd_")
            propertyName = propertyName.replace("_", " ")
            operationLongName = operationNameLookup[operationName]
            setattr(
                self,
                name,
                ct.Input(
                    doc=f"{operationLongName}-value consolidated map of {propertyName} for {coaddName} coadd",
                    name=name,
                    storageClass="HealSparseMap",
                    dimensions=("skymap", "band"),
                    multiple=False,
                    deferLoad=True,
                ),
            )


class PropertyMapSurveyWideAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=PropertyMapSurveyWideAnalysisConnections
):
    pass


class PropertyMapSurveyWideAnalysisTask(AnalysisPipelineTask):
    ConfigClass = PropertyMapSurveyWideAnalysisConfig
    _DefaultName = "propertyMapSurveyAnalysisTask"

    def parsePlotInfo(
        self, inputs: Mapping[str, Any], dataId: DataCoordinate | None, connectionNames: list[str]
    ) -> Mapping[str, Union[Mapping[str, str], str, int]]:
        """Parse the inputs and dataId to get the information needed to add to
        the figure.

        Parameters
        ----------
        inputs: `dict`
            The inputs to the task
        dataId: `~lsst.daf.butler.DataCoordinate`
            The dataId that the task is being run on.
        connectionNames: `list` [`str`]
            Name of the input connections to use for determining table names.

        Returns
        -------
        plotInfo : `dict`
            A dictionary containing the information needed to add to the
            figure.

        Notes
        -----
        We customized this method to fit our needs, because our analyses are
        not 1-1 with datasettypes. We analyze multiple connections/datasettypes
        at once, thus the table names are not the same for all connections.
        """

        # Initialize the plot info dictionary.
        plotInfo = {
            # To be filled in later.
            "tableNames": {},
            # They all share the same run, so just grab the first one.
            "run": inputs[connectionNames[0]].ref.run,
        }

        # For each connection, separately store the table name.
        for connectionName in connectionNames:
            tableName = inputs[connectionName].ref.datasetType.name
            plotInfo["tableNames"][connectionName] = tableName

        # Add the dataId information where the band, skymap and tract are the
        # same for all connections.
        self._populatePlotInfoWithDataId(plotInfo, dataId)

        return plotInfo

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        # Docstring inherited.

        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId

        data = {k: v for k, v in inputs.items() if k not in {"skymap", "healSparsePropertyMapsConfig"}}

        plotInfo = self.parsePlotInfo(inputs, dataId, list(data.keys()))
        outputs = self.run(data=data, plotConfig=self.config, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)
