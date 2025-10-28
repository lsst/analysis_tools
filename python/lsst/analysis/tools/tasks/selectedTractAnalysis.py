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

__all__ = (
    "SelectTractConfig",
    "SelectTractConnections",
)

import os

import lsst.pipe.base as pipeBase
import numpy as np
import yaml
from lsst.pex.config import Field
from lsst.pipe.base import (
    PipelineTaskConfig,
    PipelineTaskConnections,
)
from lsst.pipe.base import connectionTypes as ct


class SelectTractConnections(
    PipelineTaskConnections, dimensions=(), defaultTemplates={"metricTableName": ""}
):
    r"""Connections class to limit the scope of analysis tools tasks to the
    dataIds whose tract-level metric values lie outside the thresholds
    specified in the metricInformation.yaml file. It achieves this by using
    the adjust_all_quanta method to remove those dataIds whose metric values
    are either NaN or lie within the specified thresholds.

    This connections class should not be used standalone, instead it should
    be inherited as a second parent class. The first parent class must be an
    analysis tool connections class. The child class will first adopt the
    connections of the analysis tools connections parent class, add the
    inputTable connection defined below, then use data contained within the
    inputTable to identify the quanta that should be removed.

    For example:
    .. code-block:: python
        class SomeSelectedTractAnalysisConnections(
            SomeTractAnalysisConnections,
            SelectTractConnections,
        ):
        pass
    """

    inputTable = ct.Input(
        doc="""Table containing summary metrics for each tract.
        Used to determine which tracts exceed thresholds and
        are thus selected for plots.""",
        name="{metricTableName}",
        storageClass="ArrowAstropy",
        dimensions=("skymap",),
    )

    def adjust_all_quanta(self, adjuster):

        dataIds = list(adjuster.iter_data_ids())

        columns = ["tract"]
        if "{band}" in self.config.thresholdColumn:
            columns.extend([self.config.thresholdColumn.format(band=band) for band in self.config.bands])
            columnWithoutBand = self.config.thresholdColumn.replace("{band}_", "")
        else:
            columns.append(self.config.thresholdColumn)
            columnWithoutBand = self.config.thresholdColumn

        table = adjuster.butler.get(
            self.inputTable.name,
            dataId={"skymap": dataIds[0]["skymap"]},
            parameters={"columns": columns},
        )

        metricInfoFile = os.path.join(
            os.environ.get("ANALYSIS_TOOLS_DIR"), "metricInfo", "metricInformation.yaml"
        )
        with open(metricInfoFile) as f:
            metricInfo = yaml.load(f, Loader=yaml.SafeLoader)

        noWorkFound = True
        for dataId in dataIds:
            tractId = dataId["tract"]
            # If tract of dataId is not in table, remove quantum and move on:
            if tractId not in table["tract"]:
                adjuster.remove_quantum(dataId)
                continue

            # If metric is nan for all bands, remove quantum and move on:
            if "{band}" not in self.config.thresholdColumn:
                # Don't need to loop over bands:
                if np.isnan(table[columnWithoutBand][table["tract"] == tractId]):
                    adjuster.remove_quantum(dataId)
                    continue
            else:
                allAreNan = True
                for band in self.config.bands:
                    columnName = self.config.thresholdColumn.format(band=band)
                    if np.isfinite(table[columnName][table["tract"] == tractId]):
                        allAreNan = False
                        continue
                if allAreNan:
                    adjuster.remove_quantum(dataId)
                    continue

            # Get the relevant thresholds from the metricInfo file:
            allInThreshold = True
            columnInfo = metricInfo[columnWithoutBand]
            if "{band}" not in self.config.thresholdColumn:
                # Don't need to loop through bands:
                loThresh = columnInfo["lowThreshold"]
                hiThresh = columnInfo["highThreshold"]

                value = table[self.config.thresholdColumn][table["tract"] == tractId]
                if (value > loThresh) & (value < hiThresh):
                    adjuster.remove_quantum(dataId)
                    continue
            else:
                allInThreshold = True
                for band in self.config.bands:
                    if "lowThreshold_band" not in columnInfo:
                        loThresh = columnInfo["lowThreshold"]
                        hiThresh = columnInfo["highThreshold"]
                    else:
                        loThresh = columnInfo["lowThreshold_band"][band]
                        hiThresh = columnInfo["highThreshold_band"][band]

                    columnName = self.config.thresholdColumn.format(band=band)
                    value = table[columnName][table["tract"] == tractId]
                    if (value < loThresh) or (value > hiThresh):
                        allInThreshold = False
                        continue

            if allInThreshold:
                adjuster.remove_quantum(dataId)
                continue

            # If got this far, then there is at least one quantum to do.
            noWorkFound = False

        if noWorkFound:
            raise pipeBase.NoWorkFound(
                f"All dataIds have {self.config.thresholdColumn} values "
                "within the accepted range defined in metricInformation.yaml"
            )


class SelectTractConfig(
    PipelineTaskConfig,
    pipelineConnections=SelectTractConnections,
):
    """Config class to be used in conjunction with the SelectTractConnections
    connections class.

    This config class should not be used standalone. Instead, is should be
    inherited as a second parent class. The first parent class must be an
    analysis tool config class. The child class will first adopt config
    parameters of the analysis tools config parent class, then add to them
    the thresholdColumn config below, which is used to specify the column in
    the input metric table which should be thresholded.

    For example:
    .. code-block:: python
        class SomeSelectedTractAnalysisConfig(
            SomeTractAnalysisConfig,
            SelectTractConfig,
            pipelineConnections=SomeSelectedTractAnalysisConnections,
        ):
        pass
    """

    thresholdColumn = Field(
        doc="The name of the column in the input metric table to threshold. "
        "If it contains {band}, this will be replaced by the band name.",
        dtype=str,
        optional=False,
    )
