from __future__ import annotations

import astropy.units as apu

from typing import Iterable, cast, Mapping

from lsst.pex.config.dictField import DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField, ConfigurableActionStructField
from lsst.verify import Measurement

from ..interfaces import Tabular, TabularAction, VectorAction, MetricAction
from ..scalarActions import CountAction, MedianAction, SigmaMadAction
from ..tabularActions import AddComputedColumn, TableOfScalars, TableSelectorAction, TabularSubsetAction
from ..vectorActions import (CalcShapeSize,
                             CoaddPlotFlagSelector,
                             SnSelector,
                             StellarSelector,
                             FractionalDifference)


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


class ShapeSizeFractionalProcessBase(TabularAction):
    shapeFracDif = ConfigurableActionField(
        doc="Action which adds shapes to the table",
        default=AddComputedColumn(
            columnName="{band}_derivedShape",
            action=FractionalDifference(
                actionA=CalcShapeSize(),
                actionB=CalcShapeSize(
                    colXx="{band}_ixxPSF",
                    colYy="{band}_iyyPSF",
                    colXy="{band}_ixyPSF",
                )
            )
        ),
    )
    calculatorActions = ConfigurableActionStructField(
        doc="Actions which generate tables of scalars",
        default={
            "highSNStars": ShapeSizeFractionalScalars(
                selectors={
                    "stellar": StellarSelector(),
                    "snSelector": SnSelector(
                        threshold=2700
                    )
                }
            ),
            "lowSNStars": ShapeSizeFractionalScalars(
                selectors={
                    "stellar": StellarSelector(),
                    "snSelector": SnSelector(
                        threshold=500
                    )
                }
            )
        },
    )

    def getInputColumns(self, **kwargs) -> Iterable[str]:
        for action in self.shapeCalculator:  # type: ignore
            yield from action.getInputColumns(**kwargs)
        for action in self.calculatorActions:  # type: ignore
            yield from action.getInputColumns(**kwargs)


class ShapeSizeFractionalProcessMetric(ShapeSizeFractionalProcessBase):
    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        table = self.shapeFracDif(table, **kwargs)  # type: ignore
        results = {}
        for name, action in self.calculatorActions:  # type: ignore
            for key, value in action(table, **kwargs):
                prefix = f'{band}_' if (band := kwargs.get('band')) else ''
                newKey = f"{prefix}{name}_{key}"
                results[newKey] = value
        return results


class ShapeSizeFractionalPostProcessMetric(MetricAction):
    units = DictField(
        doc="Mapping of column name fragment to astropy unit string",
        keytype=str,
        itemtype=str,
        default={
            "highSNStars": "pixel",
            "lowSNStars": "pixel",
        }
    )

    def __call__(self, table: Tabular, **kwargs) -> Mapping[str, Measurement]:
        results = {}
        for name, value in table.items():
            # find the unit
            for fragment, unit in self.units.items():  # type:ignore
                if fragment in name:
                    break
            else:
                raise ValueError(f"Units could not be found for column {name}")
            results[name] = Measurement(name, value[0]*apu.Unit(unit))
        return results
