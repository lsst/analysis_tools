from __future__ import annotations

__all__ = (
    "DiaObjectDetectorVisitAnalysisConnections",
    "DiaObjectDetectorVisitAnalysisConfig",
    "DiaObjectDetectorVisitAnalysisTask",
)

from lsst.pipe.base import connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiaObjectDetectorVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band", "detector"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": "fakes_"},
):
    data = ct.Input(
        doc="CcdVisit-based DiaObject table to load from the butler",
        name="{fakesType}{coaddName}Diff_diaObject",
        storageClass="DataFrame",
        dimensions=("visit", "band"),
        deferLoad=True,
    )


class DiaObjectDetectorVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiaObjectDetectorVisitAnalysisConnections
):
    pass


class DiaObjectDetectorVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiaObjectDetectorVisitAnalysisConfig
    _DefaultName = "DiaObjectDetectorVisitAnalysis"
