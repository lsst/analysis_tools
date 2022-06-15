from __future__ import annotations

from typing import Optional, cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField, ConfigurableActionStructField

from ..interfaces import KeyedData, KeyedDataSchema, KeyedDataAction, Vector, VectorAction, Scalar


class KeyedDataSubsetAction(KeyedDataAction):
    columnKeys = ListField(doc="Keys to extract from KeyedData and return", dtype=str)  # type: ignore

    def getInputSchema(self) -> KeyedDataSchema:
        return ((column, Vector | Scalar) for column in self.columnKeys)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        return {key.format_map(kwargs): data[key.format_map(kwargs)] for key in self.columnKeys}  # type: ignore


class ChainedKeyedDataActions(KeyedDataAction):
    keyedDataActions = ConfigurableActionStructField(
        doc="Set of KeyedData actions to run, results will be concatenated into a final output KeyedData"
        "object"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        for action in self.tableActions:  # type: ignore
            yield from action.getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result: KeyedData = {}  # type: ignore
        for action in self.tableActions:  # type: ignore
            for column, values in action(data, **kwargs):
                result[column] = values
        return result


class AddComputedVector(KeyedDataAction):
    action = ConfigurableActionField(doc="Action to use to compute Vector", dtype=VectorAction)
    keyName = Field(doc="Key name to add to KeyedData", dtype=str)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.action.getInputSchema()  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        data[self.keyName.format(**kwargs)] = self.action(data, **kwargs)  # type: ignore
        return data


class KeyedDataSelectorAction(KeyedDataAction):
    """Down-selects rows from KeyedData of the form key: Vector

    Note this action will not work with keyed scalars, see `getInputSchema` for
    expected schema.
    """

    columnKeys = ListField(doc="Keys to extract from KeyedData and return", dtype=str)  # type: ignore

    selectors = ConfigurableActionStructField(
        doc="Selectors for selecting rows, will be AND together",
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from ((column, Vector | Scalar) for column in self.columnKeys)  # type: ignore
        for action in self.selectors:  # type: ignore
            yield from action.getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask: Optional[np.ndarray] = None
        for selector in self.selectors:  # type: ignore
            subMask = selector(data, **kwargs)
            if mask is None:
                mask = subMask
            else:
                mask *= subMask

        result = {key.format_map(kwargs): data[key.format_map(kwargs)] for key in self.columnKeys}  # type: ignore
        if mask is not None:
            return {key: cast(Vector, col)[mask] for key, col in result.items()}
        else:
            return result


class KeyedScalars(KeyedDataAction):
    scalarActions = ConfigurableActionStructField(doc="Create a KeyedData of individual ScalarActions")

    def getInputSchema(self) -> KeyedDataSchema:
        for action in self.scalarActions:  # type: ignore
            yield from action.getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result: KeyedData = {}  # type: ignore
        for name, action in self.scalarActions.items():  # type: ignore
            result[name] = action(data, **kwargs)
        return result
