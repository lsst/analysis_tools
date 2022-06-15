from __future__ import annotations

__all__ = ("BasePrep", "BaseProcess", "BaseMetricAction")

from collections import abc
from itertools import chain
from typing import cast, Mapping, Type

import astropy.units as apu

from lsst.pex.config.dictField import DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.verify import Measurement

from ..interfaces import KeyedDataAction, KeyedDataSchema, KeyedData, MetricAction, Scalar, Vector
from ..keyedDataActions import KeyedDataSelectorAction


class _PartialFormatDict(dict):
    def __missing__(self, key: str) -> str:
        return "{"+key+"}"


class BasePrep(KeyedDataSelectorAction):
    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        self.columnKeys = [name for name, _ in inputSchema]


class BaseProcess(KeyedDataAction):
    buildActions = ConfigurableActionStructField(
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
        for id, action in self.buildActions.items():  # type: ignore
            for name, typ in action.getInputSchema():
                name = name.format_map(_PartialFormatDict(identifier=id))
                buildSchema[name] = typ
        for id, action in self.filterActions.items():  # type: ignore
            for name, typ in action.getInputSchema():
                name = name.format_map(_PartialFormatDict(identifier=id))
                if name not in buildSchema:
                    filterSchema[name] = typ
        for id, action in self.calculateActions.items():  # type: ignore
            for name, typ in action.getInputSchema():
                name = name.format_map(_PartialFormatDict(identifier=id))
                if name not in buildSchema and name not in filterSchema:
                    calculateSchema[name] = typ
        return ((name, typ) for name, typ in chain(buildSchema.items(), calculateSchema.items()))

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        data = dict(data)
        for name, action in self.buildActions.items():  # type: ignore
            match action(data, **(kwargs | {'identifier': name})):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item
        view1 = data | results
        for name, action in self.filterActions.items():  # type: ignore
            match action(view1, **(kwargs | {'identifier': name})):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item

        view2 = data | results
        for name, action in self.calculateActions.items():  # type: ignore
            match action(view2, **(kwargs | {'identifier': name})):
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
