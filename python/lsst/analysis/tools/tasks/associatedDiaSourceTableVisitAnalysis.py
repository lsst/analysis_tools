from __future__ import annotations

__all__ = ("AssociatedDiaSourceTableVisitAnalysisConfig", "AssociatedDiaSourceTableVisitAnalysisTask")

from lsst.pipe.base import connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class AssociatedDiaSourceTableVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": "fakes_"},
):
    data = ct.Input(
        doc="CcdVisit-based DiaSource table to load from the butler",
        name="{fakesType}{coaddName}Diff_assocDiaSrc",
        storageClass="DataFrame",
        dimensions=("visit", "band", "detector"),
        deferLoad=True,
    )


class AssociatedDiaSourceTableVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssociatedDiaSourceTableVisitAnalysisConnections
):
    pass


class AssociatedDiaSourceTableVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AssociatedDiaSourceTableVisitAnalysisConfig
    _DefaultName = "AssociatedDiaSourceTableVisitAnalysis"
