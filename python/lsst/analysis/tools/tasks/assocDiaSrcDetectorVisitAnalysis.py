from __future__ import annotations

__all__ = ("AssocDiaSrcDetectorVisitAnalysisConfig", "AssocDiaSrcDetectorVisitAnalysisTask")

from lsst.pipe.base import connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class AssocDiaSrcDetectorVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": "fakes_"},
):
    data = ct.Input(
        doc="CcdVisit-based DiaSource table to load from the butler",
        name="{fakesType}{coaddName}Diff_assocDiaSrc",
        storageClass="DataFrame",
        dimensions=("visit", "band", "detector"),
        deferLoad=True,
    )


class AssocDiaSrcDetectorVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssocDiaSrcDetectorVisitAnalysisConnections
):
    pass


class AssocDiaSrcDetectorVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AssocDiaSrcDetectorVisitAnalysisConfig
    _DefaultName = "AssocDiaSrcDetectorVisitAnalysis"
