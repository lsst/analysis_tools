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
    "MetricAnalysisConfig",
    "MetricAnalysisTask",
)


from lsst.pex.config import ListField
from lsst.pipe.base import connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class MetricAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=(),
    defaultTemplates={"metricBundleName": ""},
):

    data = ct.Input(
        doc="A table containing metrics.",
        name="{metricBundleName}Table",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=(),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        self.dimensions.update(frozenset(sorted(config.outputDataDimensions)))
        self.data = ct.Input(
            doc=self.data.doc,
            name=self.data.name,
            storageClass=self.data.storageClass,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputDataDimensions)),
        )


class MetricAnalysisConfig(AnalysisBaseConfig, pipelineConnections=MetricAnalysisConnections):
    inputDataDimensions = ListField(
        doc="Dimensions of the input data.",
        default=(),
        dtype=str,
        optional=False,
    )
    outputDataDimensions = ListField(
        doc="Dimensions of the input data.",
        default=(),
        dtype=str,
        optional=False,
    )


class MetricAnalysisTask(AnalysisPipelineTask):
    """Perform an analysis of a metric table.
    """

    ConfigClass = MetricAnalysisConfig
    _DefaultName = "metricAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):

        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)

        data = self.loadData(inputs.pop("data"))

        if "band" in data.columns:
            outputs = self.run(data=data, plotInfo=plotInfo, bands=set(data["band"]), **inputs)
        else:
            outputs = self.run(data=data, plotInfo=plotInfo, **inputs)

        butlerQC.put(outputs, outputRefs)

