from __future__ import annotations

from typing import Mapping

from lsst.pipe.base import PipelineTask, Struct, PipelineTaskConnections, PipelineTaskConfig
from lsst.pipe.base import connectionTypes as ct
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField


from ..interfaces import KeyedData


class AnalysisBaseConnections(PipelineTaskConnections, dimensions={}):
    def __init__(self, *, config: PipelineTaskConfig = None):  # type: ignore
        super().__init__(config=config)
        existingNames = set(dir(self))
        for name in config.plots.fieldNames():  # type: ignore
            if name in self.outputs or name in existingNames:
                raise NameError(
                    f"Plot with name {name} conflicts with existing connection"
                    " are two plots named the same?"
                )
            outConnection = ct.Output(
                name=name,
                storageClass="Plot",
                doc="Dynamic connection for plotting",
                dimensions=self.dimensions,
            )
            object.__setattr__(self, name, outConnection)
            self.outputs.add(name)


class AnalysisBaseConfig(PipelineTaskConnections, connections=AnalysisBaseConnections):
    plots = ConfigurableActionStructField(doc="AnalysisPlots to run with this Task")
    metrics = ConfigurableActionStructField(doc="AnalysisMetrics to run with this Task")


class AnalysisPipelineTask(PipelineTask):
    def runPlots(self, data: KeyedData, **kwargs) -> Struct:
        results = Struct()
        for name, action in self.config.analysisPlots:  # type: ignore
            match action(data, **kwargs):
                case Mapping(val):
                    for n, v in val.items():
                        setattr(results, n, v)
                case value:
                    setattr(results, name, value)
        return results

    def runMetrics(self, data: KeyedData, **kwargs) -> Struct:
        return Struct()
