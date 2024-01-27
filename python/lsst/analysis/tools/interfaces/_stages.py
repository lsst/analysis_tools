# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
from __future__ import annotations

__all__ = ("BasePrep", "BaseProcess", "BaseMetricAction", "BaseProduce")

from collections import abc
from typing import Any, cast

import astropy.units as apu
from healsparse import HealSparseMap
from lsst.pex.config import ListField
from lsst.pex.config.configurableActions import ConfigurableActionStructField
from lsst.pex.config.dictField import DictField
from lsst.verify import Measurement

from ._actions import (
    AnalysisAction,
    JointAction,
    KeyedDataAction,
    MetricAction,
    MetricResultType,
    NoPlot,
    Tensor,
    VectorAction,
)
from ._interfaces import KeyedData, KeyedDataSchema, KeyedDataTypes, Scalar, Vector


class BasePrep(KeyedDataAction):
    """Base class for actions which prepare data for processing."""

    keysToLoad = ListField[str](doc="Keys to extract from KeyedData and return", default=[])

    vectorKeys = ListField[str](doc="Keys from the input data which selectors will be applied", default=[])

    selectors = ConfigurableActionStructField[VectorAction](
        doc="Selectors for selecting rows, will be AND together",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from (
            (column, Vector | Scalar | HealSparseMap | Tensor)
            for column in set(self.keysToLoad).union(self.vectorKeys)
        )
        for action in self.selectors:
            yield from action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return (
            (column, Vector | Scalar | HealSparseMap | Tensor)
            for column in set(self.keysToLoad).union(self.vectorKeys)
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask: Vector | None = None
        for selector in self.selectors:
            subMask = selector(data, **kwargs)
            if mask is None:
                mask = subMask
            else:
                mask *= subMask  # type: ignore
        result: dict[str, Any] = {}
        for key in set(self.keysToLoad).union(self.vectorKeys):
            formattedKey = key.format_map(kwargs)
            result[formattedKey] = cast(Vector, data[formattedKey])
        if mask is not None:
            for key in self.vectorKeys:
                # ignore type since there is not fully proper mypy support for
                # vector type casting. In the future there will be, and this
                # makes it clearer now what type things should be.
                tempFormat = key.format_map(kwargs)
                result[tempFormat] = cast(Vector, result[tempFormat])[mask]  # type: ignore
        return result

    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        existing = list(self.keysToLoad)
        existingVectors = list(self.vectorKeys)
        for name, typ in inputSchema:
            existing.append(name)
            if typ == Vector:
                existingVectors.append(name)
        self.keysToLoad = existing
        self.vectorKeys = existingVectors


class BaseProcess(KeyedDataAction):
    """Base class for actions which process data."""

    buildActions = ConfigurableActionStructField[VectorAction | KeyedDataAction](
        doc="Actions which compute a Vector which will be added to results"
    )
    filterActions = ConfigurableActionStructField[VectorAction | KeyedDataAction](
        doc="Actions which filter one or more input or build Vectors into shorter vectors"
    )
    calculateActions = ConfigurableActionStructField[AnalysisAction](
        doc="Actions which compute quantities from the input or built data"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        inputSchema: KeyedDataTypes = {}  # type: ignore
        buildOutputSchema: KeyedDataTypes = {}  # type: ignore
        filterOutputSchema: KeyedDataTypes = {}  # type: ignore
        action: AnalysisAction

        for fieldName, action in self.buildActions.items():
            for name, typ in action.getInputSchema():
                inputSchema[name] = typ
            if isinstance(action, KeyedDataAction):
                buildOutputSchema.update(action.getOutputSchema() or {})
            else:
                buildOutputSchema[fieldName] = Vector

        for fieldName, action in self.filterActions.items():
            for name, typ in action.getInputSchema():
                if name not in buildOutputSchema:
                    inputSchema[name] = typ
            if isinstance(action, KeyedDataAction):
                filterOutputSchema.update(action.getOutputSchema() or {})
            else:
                filterOutputSchema[fieldName] = Vector

        for calcAction in self.calculateActions:
            for name, typ in calcAction.getInputSchema():
                if name not in buildOutputSchema and name not in filterOutputSchema:
                    inputSchema[name] = typ
        return ((name, typ) for name, typ in inputSchema.items())

    def getOutputSchema(self) -> KeyedDataSchema:
        for action in self.buildActions:
            if isinstance(action, KeyedDataAction):
                outSchema = action.getOutputSchema()
                if outSchema is not None:
                    yield from outSchema

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        action: AnalysisAction
        results = {}
        data = dict(data)
        for name, action in self.buildActions.items():
            match action(data, **kwargs):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item
        view1 = data | results
        for name, action in self.filterActions.items():
            match action(view1, **kwargs):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item

        view2 = data | results
        for name, calcAction in self.calculateActions.items():
            match calcAction(view2, **kwargs):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item
        return results


class BaseMetricAction(MetricAction):
    """Base class for actions which compute metrics."""

    units = DictField[str, str](doc="Mapping of scalar key to astropy unit string", default={})
    newNames = DictField[str, str](
        doc="Mapping of key to new name if needed prior to creating metric",
        default={},
    )

    def getInputSchema(self) -> KeyedDataSchema:
        # Something is wrong with the typing for DictField key iteration
        return [(key, Scalar) for key in self.units]  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> MetricResultType:
        results = {}
        for key, unit in self.units.items():
            formattedKey = key.format(**kwargs)
            if formattedKey not in data:
                raise ValueError(f"Key: {formattedKey} could not be found input data")
            value = data[formattedKey]
            if not isinstance(value, Scalar):
                raise ValueError(f"Data for key {key} is not a Scalar type")
            if newName := self.newNames.get(key):
                formattedKey = newName.format(**kwargs)
            notes = {"metric_tags": kwargs.get("metric_tags", [])}
            results[formattedKey] = Measurement(formattedKey, value * apu.Unit(unit), notes=notes)
        return results


class BaseProduce(JointAction):
    """Base class for actions which produce data."""

    def setDefaults(self):
        super().setDefaults()
        self.metric = BaseMetricAction()
        self.plot = NoPlot
