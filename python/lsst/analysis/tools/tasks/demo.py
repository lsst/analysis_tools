from __future__ import annotations

from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class DemoConnections(
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


class DemoConfig(AnalysisBaseConfig, pipelineConnections=DemoConnections):
    pass


class DemoTask(AnalysisPipelineTask):
    ConfigClass = DemoConfig
    _DefaultName = "DemoAnalysis"
