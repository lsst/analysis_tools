from lsst.analysis.tools.interfaces import KeyedDataAction

from __future__ import annotations

from lsst.pex.config.listField import ListField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ..interfaces import KeyedDataAction, KeyedDataSchema, KeyedData
from ..keyedDataActions import KeyedDataSelectorAction, KeyedDataSubsetAction


class BasePrep(KeyedDataSelectorAction):
    columnKeys = ListField(doc="columns to load", dtype=str)

    def setDefaults(self):
        super().setDefaults()
        self.keyedDataAction = KeyedDataSubsetAction()
        self.keyedDataAction.columnKeys = self.columnKeys


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
