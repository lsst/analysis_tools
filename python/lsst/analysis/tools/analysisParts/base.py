from __future__ import annotations

from lsst.pex.config.listField import ListField
from lsst.pex.config.dictField import DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField
from lsst.verify import Measurement

from ..interfaces import KeyedDataAction, KeyedDataSchema, KeyedData, MetricAction, KeyedDataAction
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


class BaseMetricAction(MetricAction):
    units = DictField(
        doc="Mapping of scalar key to astropy unit string",
        keytype=str,
        itemtype=str,
    )
    newNames = DictField(
        doc="Mapping of key to new name if needed prior to creating metric",
        keytype=str,
        itemtype=str,
        default={},
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return [(cast(str, key), Scalar) for key in self.units]  # type: ignore Trouble with transitive union

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Measurement]:
        results = {}
        for key, unit in self.units.items():  # type: ignore
            formattedKey = key.format(**kwargs)
            if formattedKey not in data:
                raise ValueError(f"Key: {formattedKey} could not be found input data")
            value = data[formattedKey]
            if not isinstance(value, Scalar):
                raise ValueError(f"Data for key {key} is not a Scalar type")
            if newName := self.newNames.get(key):  # type: ignore
                formattedKey = newName.format(**kwargs)
            results[formattedKey] = Measurement(formattedKey, value * apu.Unit(unit))
        return results