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
    defaultTemplates={"metricTableName": ""},
):

    data = ct.Input(
        doc="A table containing metrics.",
        name="{metricTableName}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=(),
    )

    def __init__(self, *, config=None):

        self.dimensions.update(frozenset(sorted(config.outputDataDimensions)))
        super().__init__(config=config)
        self.data = ct.Input(
            doc=self.data.doc,
            name=self.data.name,
            storageClass=self.data.storageClass,
            deferLoad=self.data.deferLoad,
            dimensions=frozenset(sorted(config.inputDataDimensions)),
        )


class MetricAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=MetricAnalysisConnections,
):
    inputDataDimensions = ListField[str](
        doc="Dimensions of the input data table.",
        default=(),
        optional=False,
    )
    outputDataDimensions = ListField[str](
        doc="Dimensions of the outputs.",
        default=(),
        optional=False,
    )


class MetricAnalysisTask(AnalysisPipelineTask):
    """Take a metric table and run an analysis tool on the
    data it contains. This could include creating a plot
    the metrics and/or calculating summary values of those
    metrics, such as means, medians, etc. The analysis
    is outlined within the analysis tool.
    """

    ConfigClass = MetricAnalysisConfig
    _DefaultName = "metricAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Doctstring inherited

        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)

        data = self.loadData(inputs.pop("data"))

        # TODO: "bands" kwarg is a workaround for DM-47941.
        outputs = self.run(
            data=data,
            plotInfo=plotInfo,
            bands=dataId["band"],
            band=dataId["band"],
            **inputs,
        )

        butlerQC.put(outputs, outputRefs)
