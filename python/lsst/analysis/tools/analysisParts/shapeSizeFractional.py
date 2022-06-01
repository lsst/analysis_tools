from __future__ import annotations

from typing import Iterable, cast

from lsst.pipe.tasks.configurableActions import (ConfigurableActionField,
                                                 ConfigurableActionStructField)

from ..interfaces import Tabular, TabularAction, VectorAction
from ..scalarActions import CountAction, MedianAction, SigmaMadAction
from ..tabularActions import (AddComputedColumn, TableOfScalars,
                              TableSelectorAction, TabularSubsetAction)
from ..vectorActions import CalcShapeSize, CoaddPlotFlagSelector, SnSelector


class ShapeSizeFractionalPrep(TableSelectorAction):
    def setDefaults(self):
        columns = [
            "{band}_ixx",
            "{band}_iyy",
            "{band}_ixy",
            "{band}_ixxPSF",
            "{band}_iyyPSF",
            "{band}_ixyPSF",
        ]
        self.tableAction = TabularSubsetAction()
        self.tableAction.columnKeys = columns
        self.selectors.flagSelector = CoaddPlotFlagSelector()  # type: ignore
        self.selectors.snSelector = SnSelector(fluxType="psfFlux", threshold=100)  # type: ignore


class ShapeSizeFractionalScalars(TableOfScalars):
    selectors = ConfigurableActionStructField(
        doc="Selectors which determine which points go into scalar actions"
    )

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction()  # type: ignore
        self.scalarActions.sigmaMad = SigmaMadAction()  # type: ignore
        self.scalarActions.count = CountAction()  # type: ignore
        self.selectors.snSelector = SnSelector()  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        mask = kwargs.get("mask")
        for selector in cast(Iterable[VectorAction], self.selectors):  # type: ignore
            temp = selector(table, **kwargs)
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp

        return super().__call__(table, **kwargs | dict(mask=mask))


class ShapeSizeFractionalProcess(TabularAction):
    shapeCalculator = ConfigurableActionField(
        doc="Action which adds shapes to the table",
        default=AddComputedColumn(columnName="{band}_derivedShape", action=CalcShapeSize()),
    )
    calculatorActions = ConfigurableActionStructField(
        doc="Actions which generate tables of scalars",
        default={"stars": ShapeSizeFractionalScalars(), "galaxies": ShapeSizeFractionalScalars()},
    )

    def getColumns(self, **kwargs) -> Iterable[str]:
        for action in self.shapeCalculator:  # type: ignore
            yield from action.getColumns(**kwargs)
        for action in self.calculatorActions:  # type: ignore
            yield from action.getColumns(**kwargs)

    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        return super().__call__(table, **kwargs)
