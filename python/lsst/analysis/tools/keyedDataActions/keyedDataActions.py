from __future__ import annotations

__all__ = (
    "KeyedDataSubsetAction",
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

from ..interfaces import (
    KeyedData,
    KeyedDataSchema,
    KeyedDataAction,
    ScalarAction,
    Vector,
    VectorAction,
    Scalar,
)


class KeyedDataSubsetAction(KeyedDataAction):
    columnKeys = ListField[str](doc="Keys to extract from KeyedData and return")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((column, Vector | Scalar) for column in self.columnKeys)  # type: ignore

    def getOutputSchema(self) -> KeyedDataSchema:
        return ((column, Vector | Scalar) for column in self.columnKeys)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        return {key.format_map(kwargs): data[key.format_map(kwargs)] for key in self.columnKeys}  # type: ignore


class ChainedKeyedDataActions(KeyedDataAction):
    keyedDataActions = ConfigurableActionStructField[KeyedDataAction](
        doc="Set of KeyedData actions to run, results will be concatenated into a final output KeyedData"
        "object"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        for action in self.keyedDataActions:
            yield from action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        for action in self.keyedDataActions:
            yield from action.getOutputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result: KeyedData = {}  # type:ignore
        for action in self.keyedDataActions:
            for column, values in action(data, **kwargs).items():
                result[column] = values
        return result


class AddComputedVector(KeyedDataAction):
    action = ConfigurableActionField[VectorAction](doc="Action to use to compute Vector")
    keyName = Field[str](doc="Key name to add to KeyedData")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return ((self.keyName, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        data[self.keyName.format(**kwargs)] = self.action(data, **kwargs)
        return data


class KeyedDataSelectorAction(KeyedDataAction):
    """Down-selects rows from KeyedData of the form key: Vector

    Note this action will not work with keyed scalars, see `getInputSchema` for
    expected schema.
    """

    columnKeys = ListField[str](doc="Keys to extract from KeyedData and return", default=[])

    selectors = ConfigurableActionStructField[VectorAction](
        doc="Selectors for selecting rows, will be AND together",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from ((column, Vector | Scalar) for column in self.columnKeys)  # type: ignore
        for action in self.selectors:
            yield from action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        return ((column, Vector | Scalar) for column in self.columnKeys)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask: Optional[np.ndarray] = None
        for selector in self.selectors:
            subMask = selector(data, **kwargs)
            if mask is None:
                mask = subMask
            else:
                mask *= subMask  # type: ignore

        result = {key.format_map(kwargs): data[key.format_map(kwargs)] for key in self.columnKeys}  # type: ignore
        if mask is not None:
            return {key: cast(Vector, col)[mask] for key, col in result.items()}
        else:
            return result


class KeyedScalars(KeyedDataAction):
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
