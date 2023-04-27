from __future__ import annotations

__all__ = ("DiaSourceTableCcdVisitAnalysisConfig", "DiaSourceTableCcdVisitAnalysisTask")

from lsst.pipe.base import connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiaSourceTableCcdVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band"),
    defaultTemplates={"coaddName": "deep", "fakesType": "fakes_"},
):
    data = ct.Input(
        doc="CcdVisit-based DiaSource table to load from the butler",
        name="{fakesType}{coaddName}Diff_assocDiaSrc",
        storageClass="DataFrame",
        dimensions=("visit", "band", "detector"),
        deferLoad=True,
    )


class DiaSourceTableCcdVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiaSourceTableCcdVisitAnalysisConnections
):
    pass


class DiaSourceTableCcdVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiaSourceTableCcdVisitAnalysisConfig
    _DefaultName = "DiaSourceTableCcdVisitAnalysis"
