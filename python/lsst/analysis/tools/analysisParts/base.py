from lsst.analysis.tools.interfaces import KeyedDataAction

from __future__ import annotations

from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ..interfaces import KeyedDataAction, KeyedDataSchema, KeyedData


class BaseProcess(KeyedDataAction):
    metricProcess = ConfigurableActionField(doc="Does any work required to prep for metrics")
    plotProcess = ConfigurableActionField(doc="Does any work required to prep for plotting")

    @classmethod
    def getInputSchema(cls, **kwargs) -> KeyedDataSchema:
        yield from cls.plotProcess.getInputSchema(**kwargs)
        yield from cls.metricProcess.getInputSchema(**kwargs)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = dict(self.metricProcess(data, **kwargs))

        if self.plotProcess:
            results.update(self.plotProcess(dict(data) | results, **kwargs))
        return results
