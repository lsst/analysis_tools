from __future__ import annotations

from typing import Iterable, Optional, cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField, ConfigurableActionStructField

from ..interfaces import Tabular, TabularAction, Vector, VectorAction


class TabularSubsetAction(TabularAction):
    columnKeys = ListField(doc="Keys to extract from a table and return", dtype=str)  # type: ignore

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        return (column.format(**kwargs) for column in self.columnKeys)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        return {table[key] for key in self.columnKeys}  # type: ignore


class ChainedTableActions(TabularAction):
    tableActions = ConfigurableActionStructField(
        doc="Set of table actions to run, results will be concatenated into a final output table"
    )

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        for action in self.tableActions:  # type: ignore
            yield from action.getInputColumns(**kwargs)

    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        result: Tabular = {}
        for action in self.tableActions:  # type: ignore
            for column, values in action(table, **kwargs):
                result[column] = values
        return result


class AddComputedColumn(TabularAction):
    action = ConfigurableActionField(doc="Action to use to compute column", dtype=VectorAction)
    columnName = Field(doc="Column name to insert into table", dtype=str)

    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        table[self.columnName.format(**kwargs)] = self.action(table, **kwargs)  # type: ignore
        return table


class TableSelectorAction(TabularAction):
    tableAction = ConfigurableActionField(
        doc="TabularAction to run prior to down selecting rows", default=ChainedTableActions
    )
    selectors = ConfigurableActionStructField(
        doc="Selectors for selecting rows, will be AND together",
    )

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        for key in self.tableAction.getInputColumns(**kwargs):  # type: ignore
            yield key
        for action in self.selectors:  # type: ignore
            yield from action.getInputColumns(**kwargs)

    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        mask: Optional[np.ndarray] = None
        for selector in self.selectors:  # type: ignore
            subMask = selector(table, **kwargs)
            if mask is None:
                mask = subMask
            else:
                mask *= subMask

        tableResult = self.tableAction(table, **kwargs)  # type: ignore
        if mask is not None:
            return {key: col[mask] for key, col in tableResult.items()}
        else:
            return tableResult


class TableOfScalars(TabularAction):
    scalarActions = ConfigurableActionStructField(
        doc="Create a table with one row where columns correspond to individual ScalarActions"
    )

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        return (action.getInputColumns(**kwargs) for action in self.scalarActions)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        result: Tabular = {}
        for name, action in self.scalarActions.items():  # type: ignore
            result[name] = cast(Vector, np.array((action(table, **kwargs),)))
        return result
