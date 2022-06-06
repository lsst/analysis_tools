from __future__ import annotations

import astropy.units as apu
import numpy as np

from typing import Iterable, cast, Mapping, Tuple, Any

from lsst.pex.config import Field
from lsst.pex.config.dictField import DictField
from lsst.pipe.tasks.configurableActions import ConfigurableActionField, ConfigurableActionStructField
from lsst.verify import Measurement

from ..interfaces import (
    KeyedData,
    KeyedDataAction,
    Scalar,
    ScalarAction,
    VectorAction,
    MetricAction,
    KeyedDataSchema,
    Vector,
)
from ..scalarActions import CountAction, MedianAction, SigmaMadAction
from ..keyedDataActions import AddComputedVector, KeyedScalars, KeyedDataSelectorAction, KeyedDataSubsetAction
from ..vectorActions import CoaddPlotFlagSelector, SnSelector, StellarSelector, FractionalDifference
from ..vectorActions.calcShapeSize import CalcShapeSize


class _AproxMedian(ScalarAction):
    vectorKey = Field(doc="Key for the vector to perform action on", dtype=str, optional=False)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.vectorKey.format(**kwargs), Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(data[self.vectorKey.format(**kwargs)][mask])  # type: ignore
        x = int(len(value) / 10)
        return np.nanmedian(value[-x:])


class ShapeSizeFractionalScalars(KeyedScalars):
    selectors = ConfigurableActionStructField(
        doc="Selectors which determine which points go into scalar actions"
    )
    columnKey = Field(doc="Column key to compute scalars", dtype=str)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield from super().getInputSchema(**kwargs)
        for action in self.selectors:  # type: ignore
            yield from action.getInputSchema(**kwargs)

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.sigmaMad = SigmaMadAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.count = CountAction(colKey=self.columnKey)  # type: ignore
        self.selectors.snSelector = SnSelector()  # type: ignore
        self.selectors.stellar = StellarSelector()  # type: ignore
        self.scalarActions.approxMag = _AproxMedian(  # type: ignore
            vectorKey=self.selectors.snSelector.fluxType  # type: ignore
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        for selector in cast(Iterable[VectorAction], self.selectors):  # type: ignore
            temp = selector(data, **kwargs)
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp

        return super().__call__(data, **kwargs | dict(mask=mask))


class ShapeSizeFractionalPrep(KeyedDataSelectorAction):
    def setDefaults(self):
        super().setDefaults()
        # These columns must be in the output table to be used by the
        # next process step
        columns = [
            "{band}_ixx",
            "{band}_iyy",
            "{band}_ixy",
            "{band}_ixxPSF",
            "{band}_iyyPSF",
            "{band}_ixyPSF",
            "{band}_psfFlux",
            "{band}_psfFluxErr",
            "{band}_extendedness",
        ]
        self.tableAction = KeyedDataSubsetAction()
        self.tableAction.columnKeys = columns
        self.selectors.flagSelector = CoaddPlotFlagSelector()  # type: ignore
        self.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)  # type: ignore


class ShapeSizeFractionalProcessBase(KeyedDataAction):
    shapeFracDif = ConfigurableActionField(
        doc="Action which adds shapes to the KeyedData",
    )
    calculatorActions = ConfigurableActionStructField(
        doc="Actions which generate KeyedData of scalars",
    )

    def setDefaults(self):
        self.shapeFracDif = AddComputedVector()
        self.shapeFracDif.keyName = "{band}_derivedShape"
        self.shapeFracDif.action = FractionalDifference()
        self.shapeFracDif.action.actionA = CalcShapeSize()
        self.shapeFracDif.action.actionB = CalcShapeSize()
        self.shapeFracDif.action.actionB.colXx = "{band}_ixxPSF"
        self.shapeFracDif.action.actionB.colYy = "{band}_iyyPSF"
        self.shapeFracDif.action.actionB.colXy = "{band}_ixyPSF"

        self.calculatorActions.highSNStars = ShapeSizeFractionalScalars(
            columnKey=self.shapeFracDif.columnName
        )
        self.calculatorActions.lowSNStars = ShapeSizeFractionalScalars(columnKey=self.shapeFracDif.columnName)

        for action in self.calculatorActions:  # type: ignore
            action.setDefaults()

        self.calculatorActions.highSNStars.selectors
        self.calculatorActions.highSNStars.selectors.stellar = SnSelector(threshold=2700)

        self.calculatorActions.lowSNStars.selectors.update = {  # type: ignore
            "stellar": StellarSelector(),
            "snSelector": SnSelector(threshold=500),
        }
        super().setDefaults()

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield from self.shapeFracDif.getInputSchema(**kwargs)  # type: ignore
        for action in self.calculatorActions:  # type: ignore
            yield from action.getInputSchema(**kwargs)


class ShapeSizeFractionalProcessMetric(ShapeSizeFractionalProcessBase):
    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        table = self.shapeFracDif(data, **kwargs)  # type: ignore
        results = {}
        for name, action in self.calculatorActions.items():  # type: ignore
            for key, value in action(table, **kwargs).items():
                prefix = f"{band}_" if (band := kwargs.get("band")) else ""
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
        },
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return [(cast(str, key), Scalar) for key in self.units]  # type: ignore Trouble with transitive union

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Measurement]:
        results = {}
        for name, value in data.items():
            # find the unit
            for fragment, unit in self.units.items():  # type:ignore
                if fragment in name:
                    break
            else:
                raise ValueError(f"Units could not be found for column {name}")
            results[name] = Measurement(name, value * apu.Unit(unit))
        return results
