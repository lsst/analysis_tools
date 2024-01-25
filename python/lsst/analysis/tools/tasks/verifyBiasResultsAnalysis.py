from __future__ import annotations

__all__ = (
    "VerifyBiasResultsAnalysisConnections",
    "VerifyBiasResultsAnalysisConfig",
    "VerifyBiasResultsAnalysisTask",
)

from lsst.pipe.base.connectionTypes import Input

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class VerifyBiasResultsAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument",),
):
    data = Input(
        doc="verifyBiasResults",
        name="verifyBiasResults",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
    )


class VerifyBiasResultsAnalysisConfig(
    AnalysisBaseConfig,
    pipelineConnections=VerifyBiasResultsAnalysisConnections,
):
    pass


class VerifyBiasResultsAnalysisTask(AnalysisPipelineTask):
    ConfigClass = VerifyBiasResultsAnalysisConfig
    _DefaultName = "verifyBiasResultsAnalysis"
