from __future__ import annotations

from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConnections, AnalysisBaseConfig, AnalysisPipelineTask

from ..analysisPlots.analysisPlots import ShapeSizeFractionalDiffScatter
from ..analysisMetrics.analysisMetrics import ShapeSizeFractionalMetric

from ..vectorActions.selectors import VisitPlotFlagSelector

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
    )


class SourceTableVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=SourceTableVisitAnalysisConnections
):

    def setDefaults(self):
        super().setDefaults()
        # set plots to run
        self.plots.shapeSizeFractionalDiffScatter = ShapeSizeFractionalDiffScatter()
        self.plots.shapeSizeFractionalDiffScatter.flagSelector = VisitPlotFlagSelector()

        # set metrics to run
        self.metrics.shapeSizeFractionalMetric = ShapeSizeFractionalMetric()


class sourceTableVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SourceTableVisitAnalysisConfig
    _DefaultName = "sourceTableVisitAnalysisTask"
