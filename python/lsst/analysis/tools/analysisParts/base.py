from __future__ import annotations

__all__ = ("BasePrep", "BaseProcess", "BaseMetricAction")

from collections import abc
from itertools import chain
from typing import cast, Mapping, Type

import astropy.units as apu

from lsst.pex.config.dictField import DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.verify import Measurement

from ..interfaces import AnalysisAction, KeyedDataAction, KeyedDataSchema, KeyedData, MetricAction, Scalar, Vector
from ..keyedDataActions import KeyedDataSelectorAction


class BasePrep(KeyedDataSelectorAction):
    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        self.columnKeys = [name for name, _ in inputSchema]


class BaseProcess(KeyedDataAction):
    buildActions = ConfigurableActionStructField[AnalysisAction](
        doc="Actions which compute a Vector which will be added to results"
    )
    filterActions = ConfigurableActionStructField(
        doc="Actions which filter one or more input or build Vectors into shorter vectors"
    )
    calculateActions = ConfigurableActionStructField(
        doc="Actions which compute quantities from the input or built data"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        buildSchema: dict[str, Type[Vector] | Type[Scalar]] = {}
        filterSchema: dict[str, Type[Vector] | Type[Scalar]] = {}
        calculateSchema: dict[str, Type[Vector] | Type[Scalar]] = {}
        for action in self.buildActions:
            for name, typ in action.getInputSchema():
                buildSchema[name] = typ
        for action in self.filterActions:  # type: ignore
            for name, typ in action.getInputSchema():
                if name not in buildSchema:
                    filterSchema[name] = typ
        for action in self.calculateActions:  # type: ignore
            for name, typ in action.getInputSchema():
                if name not in buildSchema and name not in filterSchema:
                    calculateSchema[name] = typ
        return ((name, typ) for name, typ in chain(buildSchema.items(), calculateSchema.items()))

    def getOutputSchema(self) -> KeyedDataSchema:
        for action in self.buildActions:
            if isinstance(action, KeyedDataAction):
                yield from action.getOutputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        data = dict(data)
        for name, action in self.buildActions.items():  # type: ignore
            match action(data, **kwargs):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item
        view1 = data | results
        for name, action in self.filterActions.items():  # type: ignore
            match action(view1, **kwargs):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item

        view2 = data | results
        for name, action in self.calculateActions.items():  # type: ignore
            match action(view2, **kwargs):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item
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

    def getInputSchema(self) -> KeyedDataSchema:
        return [(cast(str, key), Scalar) for key in self.units]  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Measurement] | Measurement:
        results = {}
        for key, unit in self.units.items():  # type: ignore
            formattedKey = key.format(**kwargs)
            if formattedKey not in data:
                raise ValueError(f"Key: {formattedKey} could not be found input data")
            value = data[formattedKey]
            if not isinstance(value, Scalar):
                import ipdb;ipdb.set_trace()
                raise ValueError(f"Data for key {key} is not a Scalar type")
            if newName := self.newNames.get(key):  # type: ignore
                formattedKey = newName.format(**kwargs)
            results[formattedKey] = Measurement(formattedKey, value * apu.Unit(unit))
        return results
