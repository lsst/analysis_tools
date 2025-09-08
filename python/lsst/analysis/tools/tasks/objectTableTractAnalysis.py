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
    "ObjectTableTractAnalysisConnections",
    "ObjectTableTractAnalysisConfig",
    "ObjectTableTractAnalysisTask",
)

import copy
from typing import TYPE_CHECKING

import lsst.pex.config as pexConfig
from lsst.obs.base.utils import TableVStack
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

if TYPE_CHECKING:
    # from lsst.daf.butler import DatasetRef, DeferredDatasetHandle
    from lsst.pipe.base.connections import InputQuantizedConnection, OutputQuantizedConnection
    from lsst.pipe.base import QuantumContext

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class ObjectTableTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={"inputName": "objectTable_tract", "outputName": "objectTable_tract"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="{inputName}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("skymap", "tract"),
    )

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    def __init__(self, *, config=None):
        if config:
            if not config.load_skymap:
                del self.skymap
            if config.load_multiple:
                connection = self.data
                del self.data
                self.data = type(connection)(
                    doc=connection.doc,
                    name=connection.name,
                    storageClass=connection.storageClass,
                    deferLoad=connection.deferLoad,
                    dimensions=connection.dimensions,
                    multiple=True,
                )
                self.dimensions = {"skymap"}
        super().__init__(config=config)


class ObjectTableTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=ObjectTableTractAnalysisConnections
):
    load_multiple = pexConfig.Field[bool](
        doc="Whether to load multiple tables and concatenate them to make one output for all tracts.",
        default=False,
    )
    load_skymap = pexConfig.Field[bool](doc="Whether to load the skymap", default=True)


class ObjectTableTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = ObjectTableTractAnalysisConfig
    _DefaultName = "objectTableTractAnalysis"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        if not self.config.load_multiple:
            super().runQuantum(butlerQC, inputRefs, outputRefs)
            return
        inputs = butlerQC.get(inputRefs)
        extra_values = {}
        handles = {}
        columns = self.collectInputNames()
        add_tract = "tract" not in columns

        for tract, handle, idx_h in sorted(
            (inputRef.dataId["tract"], inputHandle, idx)
            for idx, (inputRef, inputHandle) in enumerate(zip(inputRefs.data, inputs["data"], strict=True))
        ):
            handles[tract] = handle
            if add_tract:
                extra_values[idx_h] = {"tract": tract}
        data = TableVStack.vstack_handles(
            handles.values(), extra_values=extra_values, kwargs_get={"parameters": {"columns": columns}}
        )

        tracts = ",".join(str(tract) for tract in sorted(handles.keys()))
        input_first = inputs.pop("data")[0]
        inputs_plot = {"data": input_first}
        dataId = copy.copy(input_first.dataId)
        dataId._values = (dataId._values[0], tracts)
        plotInfo = self.parsePlotInfo(inputs_plot, dataId)

        outputs = self.run(data=data, plotInfo=plotInfo, **inputs)
        self.putByBand(butlerQC, outputs, outputRefs)
