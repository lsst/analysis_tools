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

__all__ = ("MetadataAnalysisConfig", "MetadataAnalysisTask")

import lsst.pex.config
from lsst.pipe.base import NoWorkFound, connectionTypes

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class MetadataAnalysisConnections(
    AnalysisBaseConnections,
    dimensions={"instrument", "exposure", "detector"},
    defaultTemplates={"taskName": "isr"},
):
    data = connectionTypes.Input(
        doc="Task metadata to load.",
        name="{taskName}_metadata",
        storageClass="TaskMetadata",
        deferLoad=True,
        dimensions={"instrument", "exposure", "detector"},
    )

    def __init__(self, *, config=None):
        """Customize the connections for a specific task's metadata.

        Parameters
        ----------
        config : `MetadataMetricConfig`
            A config for `MetadataMetricTask` or one of its subclasses.
        """
        super().__init__(config=config)
        if config and config.metadataDimensions != self.data.dimensions:
            self.dimensions.clear()
            self.dimensions.update(config.metadataDimensions)
            self.data = connectionTypes.Input(
                name=self.data.name,
                doc=self.data.doc,
                storageClass=self.data.storageClass,
                deferLoad=self.data.deferLoad,
                dimensions=frozenset(config.metadataDimensions),
            )


class MetadataAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=MetadataAnalysisConnections,
):
    metadataDimensions = lsst.pex.config.ListField(
        # Sort to ensure default order is consistent between runs
        default=sorted(MetadataAnalysisConnections.dimensions),
        dtype=str,
        doc="Override for the dimensions of the input data.",
    )


class MetadataAnalysisTask(AnalysisPipelineTask):
    ConfigClass = MetadataAnalysisConfig
    _DefaultName = "metadataAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        dataId = butlerQC.quantum.dataId
        inputs = butlerQC.get(inputRefs)
        plotInfo = self.parsePlotInfo(inputs, dataId)
        metadata = inputs["data"].get().to_dict()
        if not metadata:
            taskName = inputRefs.data.datasetType.name
            taskName = taskName[: taskName.find("_")]
            raise NoWorkFound(f"No metadata entries for {taskName}.")
        outputs = self.run(data=metadata, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)
