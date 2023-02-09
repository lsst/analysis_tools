from __future__ import annotations

__all__ = (
    "DiaObjectTableAnalysisConnections",
    "DiaObjectTableAnalysisConfig",
    "DiaObjectTableAssociatedSourcesTask",
)

from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DiaObjectTableAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band"),
    defaultTemplates={"coaddName": "deep", "fakesType": ""},
):
    data = ct.Input(
        doc="CcdVisit-based DiaObject table to load from the butler",
        name="{fakesType}{coaddName}Diff_diaObject",
        storageClass="DataFrame",
        dimensions=("visit", "band", "detector"),
        deferLoad=True,
    )


class DiaObjectTableAnalysisConfig(AnalysisBaseConfig, pipelineConnections=DiaObjectTableAnalysisConnections):
    pass


class DiaObjectTableAssociatedSourcesTask(AnalysisPipelineTask):
    ConfigClass = DiaObjectTableAnalysisConfig
    _DefaultName = "DiaObjectTableAssociatedSources"
