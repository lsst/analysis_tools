from __future__ import annotations

from lsst.pipe.base import connectionTypes as ct

from ..analysisMetrics.analysisMetrics import ShapeSizeFractionalMetric
from ..analysisPlots.analysisPlots import ShapeSizeFractionalDiffScatter

# from ..vectorActions.selectors import VisitPlotFlagSelector
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
    )


class SourceTableVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=SourceTableVisitAnalysisConnections
):
    def setDefaults(self):
        super().setDefaults()
        # set plots to run
        self.plots.shapeSizeFractionalDiffScatter = ShapeSizeFractionalDiffScatter()
        # self.plots.shapeSizeFractionalDiffScatter.flagSelector = \
        # VisitPlotFlagSelector()

        # set metrics to run
        self.metrics.shapeSizeFractionalMetric = ShapeSizeFractionalMetric()


class sourceTableVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SourceTableVisitAnalysisConfig
    _DefaultName = "sourceTableVisitAnalysisTask"
