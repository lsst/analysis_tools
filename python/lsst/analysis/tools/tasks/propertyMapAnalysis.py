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
    "SurveyWidePropertyMapAnalysisConfig",
    "SurveyWidePropertyMapAnalysisTask",
    "PerTractPropertyMapAnalysisConfig",
    "PerTractPropertyMapAnalysisTask",
]

from typing import Any, Mapping, Union

from lsst.daf.butler import DataCoordinate
from lsst.pex.config import ChoiceField, DictField, Field, ListField
from lsst.pipe.base import (
    InputQuantizedConnection,
    OutputQuantizedConnection,
    QuantumContext,
    connectionTypes,
)
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class PerTractPropertyMapAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "band", "tract"),
    defaultTemplates={"outputName": "propertyMapTract"},
):
    skymap = connectionTypes.Input(
        doc="The skymap that covers the tract that the data is from.",
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
            propertyName, operationName = name.split("_map_")
            coaddName, propertyName = propertyName.split("Coadd_")
            propertyName = propertyName.replace("_", " ")
            operationLongName = operationNameLookup[operationName]
            setattr(
                self,
                name,
                connectionTypes.Input(
                    doc=f"{operationLongName}-value map of {propertyName} for {coaddName} coadd",
                    name=name,
                    storageClass="HealSparseMap",
                    dimensions=("skymap", "band", "tract"),
                    multiple=False,
                    deferLoad=True,
                ),
            )


class PerTractPropertyMapAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=PerTractPropertyMapAnalysisConnections
):
    projectionKwargs = DictField(
        keytype=str,
        itemtype=float,
        doc="Keyword arguments to use in the GnomonicSkyproj call, e.g. n_grid_lon. "
        "The following keys are not permitted: 'ax', 'lon_0', 'lat_0', and 'extent'. "
        "See https://skyproj.readthedocs.io/en/latest/modules.html#skyproj.skyproj.GnomonicSkyproj",
        default={},
        keyCheck=lambda k: k not in ["ax", "lon_0", "lat_0", "extent"],
    )

    zoomFactors = ListField(
        dtype=float,
        doc="Two-element list of zoom factors to use when plotting the maps.",
        default=[2, 8],
    )

    colorbarKwargs = DictField(
        keytype=str,
        itemtype=str,
        doc="Keyword arguments to pass to the colorbar except for 'orientation' and 'location'.",
        default={"cmap": "viridis"},
        keyCheck=lambda k: k not in ["orientation", "location"],
    )

    publicationStyle = Field[bool](doc="Make a simplified plot for publication use?", default=False)


class PerTractPropertyMapAnalysisTask(AnalysisPipelineTask):
    ConfigClass = PerTractPropertyMapAnalysisConfig
    _DefaultName = "perTractPropertyMapAnalysisTask"

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

        # Add the dataId information, same for all connections.
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
        skymap = inputs["skymap"]
        dataId = butlerQC.quantum.dataId
        tractInfo = skymap[dataId["tract"]]
        mapKeys = [key for key in inputs if key != "skymap"]

        plotInfo = self.parsePlotInfo(inputs, dataId, mapKeys)
        outputs = self.run(data=inputs, tractInfo=tractInfo, plotConfig=self.config, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)


class SurveyWidePropertyMapAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "band"),
    defaultTemplates={"outputName": "propertyMapSurvey"},
):
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
                connectionTypes.Input(
                    doc=f"{operationLongName}-value consolidated map of {propertyName} for {coaddName} coadd",
                    name=name,
                    storageClass="HealSparseMap",
                    dimensions=("skymap", "band"),
                    multiple=False,
                    deferLoad=True,
                ),
            )


class SurveyWidePropertyMapAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=SurveyWidePropertyMapAnalysisConnections
):
    # Note: Gnomonic projection is excluded here because `GnomonicSkyproj` must
    # have the central lon/lat set (defaults to 0/0) which makes it useful for
    # plotting individual tracts but not for survey-wide maps.
    projection = ChoiceField[str](
        doc="The projection to use for plotting the map. "
        "See https://skyproj.readthedocs.io/en/latest/projections.html",
        default="McBryde",
        allowed={
            proj: proj
            for proj in (
                "McBryde",
                "Mollweide",
                "Cylindrical",
                "Laea",
                "Hammer",
                "EqualEarth",
                "ObliqueMollweide",
                "Albers",
            )
        },
    )

    projectionKwargs = DictField(
        keytype=str,
        itemtype=float,
        doc="Keyword arguments to use in the projection call, e.g. lon_0. The key 'ax' is not permitted. "
        "See https://skyproj.readthedocs.io/en/latest/projections.html",
        default={},
        keyCheck=lambda k: k not in ["ax"],
    )

    autozoom = Field[bool](
        doc="Automatically zooms in on the RA/Dec range of the map to make better use of its resolution; "
        "otherwise, the map is displayed within the full-sky domain.",
        default=True,
    )

    colorbarKwargs = DictField(
        keytype=str,
        itemtype=str,
        doc="Keyword arguments to pass to the colorbar.",
        default={"orientation": "horizontal", "location": "top", "cmap": "viridis"},
    )
    # TODO: Mixed types will be allowed after DM-47937. You can then add things
    # like "aspect": 20, rather than just strings.

    publicationStyle = Field[bool](doc="Make a simplified plot for publication use?", default=False)


class SurveyWidePropertyMapAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SurveyWidePropertyMapAnalysisConfig
    _DefaultName = "surveyWidePropertyMapAnalysisTask"

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

        # Add the dataId information, same for all connections.
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

        plotInfo = self.parsePlotInfo(inputs, dataId, list(inputs.keys()))
        outputs = self.run(data=inputs, plotConfig=self.config, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)
