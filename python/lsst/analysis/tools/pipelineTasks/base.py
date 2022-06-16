from __future__ import annotations

from collections import abc
from typing import cast, Iterable

from lsst.pipe.base import PipelineTask, Struct, PipelineTaskConnections, PipelineTaskConfig
from lsst.pipe.base import connectionTypes as ct
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.pex.config import ListField


from ..interfaces import KeyedData
from ..analysisMetrics.metricActionMapping import MetricActionMapping


class AnalysisBaseConnections(
        PipelineTaskConnections,
        dimensions={},
        defalutTemplates={"inputName": "Placeholder"}
        ):

    metrics = ct.Output(
        doc="Metrics calculated on input dataset type",
        name="{inputName}_metrics",
        storageClass=""
    )


    def __init__(self, *, config: PipelineTaskConfig = None):  # type: ignore
        if (inputName := config.connections.inputName) == "Placeholder":  # type: ignore
            raise RuntimeError("Subclasses must specify an alternative value for the defaultTemplate `inputName`")
        super().__init__(config=config)
        existingNames = set(dir(self))
        for name in config.plots.fieldNames():  # type: ignore
            name = f"{inputName}_{name}"
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

    bands = ListField(
        doc="Filter bands on which to run all of the actions",
        dtype=str,
        default=["g", "r", "i", "z", "y"]
    )


class AnalysisPipelineTask(PipelineTask):
    def runPlots(self, data: KeyedData, **kwargs) -> Struct:
        results = Struct()
        for name, action in self.config.analysisPlots.items():  # type: ignore
            match action(data, **kwargs):
                case abc.Mapping() as val:
                    for n, v in val.items():
                        setattr(results, n, v)
                case value:
                    setattr(results, name, value)
        return results

    def runMetrics(self, data: KeyedData, **kwargs) -> Struct:
        metricsMapping = MetricActionMapping()
        for name, action in self.config.metrics.items():  # type: ignore
            match action(data, **kwargs):
                case abc.Mapping() as val:
                    results = list(val.values())
                case val:
                    results = [val]
            metricsMapping[name] = results
        return Struct(metrics=metricsMapping) 

    def run(self, data: KeyedData, **kwargs) -> Struct:
        results = Struct()

        kwargs['bands'] = cast(Iterable[str], self.config.bands):  # type: ignore
        results.mergeItems(self.runPlots(data, **kwargs))
        results.mergeItems(self.runMetrics(data, **kwargs))

        return results
