from __future__ import annotations

from collections import abc
from lsst.pipe.base import PipelineTask, Struct, PipelineTaskConnections, PipelineTaskConfig
from lsst.pipe.base import connectionTypes as ct
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.pex.config import ListField


from ..interfaces import KeyedData, AnalysisPlot, AnalysisMetric
from ..analysisMetrics.metricActionMapping import MetricActionMapping


class AnalysisBaseConnections(
    PipelineTaskConnections, dimensions={}, defaultTemplates={"inputName": "Placeholder"}
):

    metrics = ct.Output(
        doc="Metrics calculated on input dataset type",
        name="{inputName}_metrics",
        storageClass="AnalysisMetricStack"
    )

    def __init__(self, *, config: "AnalysisBaseConfig" = None):  # type: ignore
        if (inputName := config.connections.inputName) == "Placeholder":  # type: ignore
            raise RuntimeError(
                "Subclasses must specify an alternative value for the defaultTemplate `inputName`"
            )
        super().__init__(config=config)

        # Set the dimensions for the metric
        self.metrics = ct.Output(
            name=self.metrics.name,
            doc=self.metrics.doc,
            storageClass=self.metrics.storageClass,
            dimensions=self.dimensions,
            multiple=False,
            isCalibration=False
        )

        existingNames = set(dir(self))
        names = []
        for plotAction in config.plots:
            if plotAction.multiband:
                for band in config.bands:
                    names.extend(name.format(band=band) for name in plotAction.getOutputNames())
            else:
                names.extend(plotAction.getOutputNames())
        for name in names:
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


class AnalysisBaseConfig(PipelineTaskConfig, pipelineConnections=AnalysisBaseConnections):
    plots = ConfigurableActionStructField[AnalysisPlot](doc="AnalysisPlots to run with this Task")
    metrics = ConfigurableActionStructField[AnalysisMetric](doc="AnalysisMetrics to run with this Task")

    bands = ListField[str](
        doc="Filter bands on which to run all of the actions", default=["g", "r", "i", "z", "y"]
    )


class StandinPlotInfo(dict):
    def __missing__(self, key):
        return ""


class AnalysisPipelineTask(PipelineTask):
    # Typing config because type checkers dont know about our Task magic
    config: AnalysisBaseConfig
    ConfigClass = AnalysisBaseConfig

    def runPlots(self, data: KeyedData, **kwargs) -> Struct:
        results = Struct()
        # allow not sending in plot info
        if 'plotInfo' not in kwargs:
            kwargs['plotInfo'] = StandinPlotInfo()
        for name, action in self.config.plots.items():
            match action(data, **kwargs):
                case abc.Mapping() as val:
                    for n, v in val.items():
                        setattr(results, n, v)
                case value:
                    setattr(results, name, value)
        return results

    def runMetrics(self, data: KeyedData, **kwargs) -> Struct:
        metricsMapping = MetricActionMapping()
        for name, action in self.config.metrics.items():
            match action(data, **kwargs):
                case abc.Mapping() as val:
                    results = list(val.values())
                case val:
                    results = [val]
            metricsMapping[name] = results  # type: ignore
        return Struct(metrics=metricsMapping)

    def run(self, data: KeyedData, **kwargs) -> Struct:
        results = Struct()
        plotKey = f"{self.config.connections.inputName}_{{name}}"  # type: ignore
        if "bands" not in kwargs:
            kwargs["bands"] = list(self.config.bands)
        for name, value in self.runPlots(data, **kwargs).getDict().items():
            setattr(results, plotKey.format(name=name), value)
        for name, value in self.runMetrics(data, **kwargs).getDict().items():
            setattr(results, name, value)

        return results
