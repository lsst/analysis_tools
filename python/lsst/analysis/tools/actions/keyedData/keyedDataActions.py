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

__all__ = (
    "ChainedKeyedDataActions",
    "AddComputedVector",
    "KeyedDataSelectorAction",
    "KeyedScalars",
)

from typing import Optional, cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField, ConfigurableActionStructField

from ...interfaces import (
    KeyedData,
    KeyedDataAction,
    KeyedDataSchema,
    Scalar,
    ScalarAction,
    Vector,
    VectorAction,
)


class ChainedKeyedDataActions(KeyedDataAction):
    r"""Run a series of `KeyedDataAction`\ s and accumulated their output into
    one KeyedData result.
    """

    keyedDataActions = ConfigurableActionStructField[KeyedDataAction](
        doc="Set of KeyedData actions to run, results will be concatenated into a final output KeyedData"
        "object"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        for action in self.keyedDataActions:
            yield from action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        for action in self.keyedDataActions:
            output = action.getOutputSchema()
            if output is not None:
                yield from output

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result: KeyedData = {}  # type:ignore
        for action in self.keyedDataActions:
            for column, values in action(data, **kwargs).items():
                result[column] = values
        return result


class AddComputedVector(KeyedDataAction):
    """Compute a `Vector` from the specified `VectorAction` and add it to a
    copy of the KeyedData, returning the result.
    """

    action = ConfigurableActionField[VectorAction](doc="Action to use to compute Vector")
    keyName = Field[str](doc="Key name to add to KeyedData")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return ((self.keyName, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result = dict(data)
        result[self.keyName.format(**kwargs)] = self.action(data, **kwargs)
        return result


class KeyedDataSelectorAction(KeyedDataAction):
    """Extract Vector specified by ``vectorKeys`` from input KeyedData and
    optionally apply selectors to down select extracted vectors.

    Note this action will not work with keyed scalars, see `getInputSchema` for
    expected schema.
    """

    vectorKeys = ListField[str](doc="Keys to extract from KeyedData and return", default=[])

    selectors = ConfigurableActionStructField[VectorAction](
        doc="Selectors for selecting rows, will be AND together",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from ((column, Vector | Scalar) for column in self.vectorKeys)  # type: ignore
        for action in self.selectors:
            yield from action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return ((column, Vector | Scalar) for column in self.vectorKeys)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask: Optional[np.ndarray] = None
        for selector in self.selectors:
            subMask = selector(data, **kwargs)
            if mask is None:
                mask = subMask
            else:
                mask *= subMask  # type: ignore

        result = {key.format_map(kwargs): data[key.format_map(kwargs)] for key in self.vectorKeys}
        if mask is not None:
            return {key: cast(Vector, col)[mask] for key, col in result.items()}
        else:
            return result


class KeyedScalars(KeyedDataAction):
    """Creates an output of type KeyedData, where the keys are given byt the
    identifiers in `scalarActions` and the values are the results of the
    corresponding `ScalarAction`.
    """

    scalarActions = ConfigurableActionStructField[ScalarAction](
        doc="Create a KeyedData of individual ScalarActions"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        for action in self.scalarActions:
            yield from action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        for name in self.scalarActions.fieldNames:
            yield (name, Scalar)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result: KeyedData = {}  # type: ignore
        for name, action in self.scalarActions.items():
            result[name] = action(data, **kwargs)
        return result
