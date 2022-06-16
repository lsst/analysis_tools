from __future__ import annotations


from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConnections, AnalysisBaseConfig, AnalysisPipelineTask


class ObjectTableTractAnalysisConnections(AnalysisBaseConnections):
    objectTable = ct.Input(
        doc="Tract based object table to load from the butler",
        name="objectTable_tract",
        storageClass="DataFrame",
        deferLoad=True
    )
    metrics = ct.Output(
        doc="Metrics calculated on tract based object table",
        name="objectTableTract_metrics"
    )


class ObjectTableTractAnalysisConfig(AnalysisBaseConfig):
    pass