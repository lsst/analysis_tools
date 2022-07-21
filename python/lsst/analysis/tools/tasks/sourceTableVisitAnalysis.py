from __future__ import annotations

from lsst.pipe.base import connectionTypes as ct

from ..analysisMetrics import SkyFluxStatisticMetric
from ..analysisMetrics.limitingMagnitudeMetric import FiveSigmaPointSourceDepthMetric
from ..contexts import VisitContext
from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class SourceTableVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band"),
    defaultTemplates={"inputName": "sourceTable_visit"},
):

    data = ct.Input(
        doc="Visit based source table to load from the butler",
        name="sourceTable_visit",
        storageClass="DataFrame",
        dimensions=("visit", "band"),
        deferLoad=True,
    )


class SourceTableVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=SourceTableVisitAnalysisConnections
):
    def setDefaults(self):
        super().setDefaults()
        # set plots to run

        # set metrics to run
        self.metrics.fiveSigmaPointSourceDepthMetric = FiveSigmaPointSourceDepthMetric()
        self.metrics.skyFluxVisitStatisticMetric = SkyFluxStatisticMetric()
        self.metrics.skyFluxVisitStatisticMetric.applyContext(VisitContext)


class SourceTableVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SourceTableVisitAnalysisConfig
    _DefaultName = "sourceTableVisitAnalysisTask"
