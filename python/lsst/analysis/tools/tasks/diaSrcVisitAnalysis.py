from __future__ import annotations

__all__ = ("DiaSrcCcdVisitAnalysisConfig", "DiaSrcCcdVisitAnalysisTask")

from lsst.pipe.base import connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiaSrcCcdVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band"),
    defaultTemplates={"coaddName": "deep"},
):
    data = ct.Input(
        doc="CcdVisit-based DiaSource table to load from the butler",
        name="{coaddName}Diff_diaSrc",
        storageClass="SourceCatalog",
        dimensions=("visit", "band", "detector"),
        deferLoad=True,
    )


class DiaSrcCcdVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiaSrcCcdVisitAnalysisConnections
):
    pass


class DiaSrcCcdVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiaSrcCcdVisitAnalysisConfig
    _DefaultName = "DiaSrcCcdVisitAnalysis"
