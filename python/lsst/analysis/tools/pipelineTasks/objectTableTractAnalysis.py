from __future__ import annotations


from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConnections, AnalysisBaseConfig, AnalysisPipelineTask

from ..analysisPlots.analysisPlots import ShapeSizeFractionalDiffScatter
from ..analysisMetrics.analysisMetrics import ShapeSizeFractionalMetric


class ObjectTableTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={"inputName": "objectTable_tract"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="objectTable_tract",
        storageClass="DataFrame",
        #deferLoad=True,
        dimensions=("skymap", "tract"),
    )


class ObjectTableTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=ObjectTableTractAnalysisConnections
):
    def setDefaults(self):
        super().setDefaults()
        # set plots to run
        self.plots.shapeSizeFractionalDiffScatter = ShapeSizeFractionalDiffScatter()

        # set metrics to run
        self.metrics.shapeSizeFractionalMetric = ShapeSizeFractionalMetric()


class ObjectTableTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = ObjectTableTractAnalysisConfig
    _DefaultName = "objectTableTractAnalysisTask"
