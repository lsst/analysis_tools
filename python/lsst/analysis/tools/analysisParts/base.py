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

__all__ = ("BasePrep", "BaseProcess", "BaseMetricAction")

from collections import abc
from typing import Mapping

import astropy.units as apu
from lsst.pex.config.dictField import DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionStructField
from lsst.verify import Measurement

from ..actions.keyedData import KeyedDataSelectorAction
from ..interfaces import (
    AnalysisAction,
    KeyedData,
    KeyedDataAction,
    KeyedDataSchema,
    KeyedDataTypes,
    MetricAction,
    Scalar,
    Vector,
    VectorAction,
)


class BasePrep(KeyedDataSelectorAction):
    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        self._frozen = False
        self.vectorKeys = [name for name, _ in inputSchema]
        self._frozen = True


class BaseProcess(KeyedDataAction):
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

        for action in self.calculateActions:
            for name, typ in action.getInputSchema():
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
        for name, action in self.calculateActions.items():
            match action(view2, **kwargs):
                case abc.Mapping() as item:
                    for key, result in item.items():
                        results[key] = result
                case item:
                    results[name] = item
        return results


class BaseMetricAction(MetricAction):
    units = DictField[str, str](
        doc="Mapping of scalar key to astropy unit string",
    )
    newNames = DictField[str, str](
        doc="Mapping of key to new name if needed prior to creating metric",
        default={},
    )

    def getInputSchema(self) -> KeyedDataSchema:
        # Something is wrong with the typing for DictField key iteration
        return [(key, Scalar) for key in self.units]  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Measurement] | Measurement:
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
            results[formattedKey] = Measurement(formattedKey, value * apu.Unit(unit))
        return results
