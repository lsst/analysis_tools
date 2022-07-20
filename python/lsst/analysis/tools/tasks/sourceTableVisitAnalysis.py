from __future__ import annotations

from lsst.pipe.base import connectionTypes as ct

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
    pass


class SourceTableVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SourceTableVisitAnalysisConfig
    _DefaultName = "sourceTableVisitAnalysis"
