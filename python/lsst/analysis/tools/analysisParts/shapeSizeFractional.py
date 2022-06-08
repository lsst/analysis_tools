from __future__ import annotations

import astropy.units as apu
import numpy as np

from typing import cast, Mapping
from lsst.analysis.tools.plotActions.scatterplotWithTwoHists import ScatterPlotWithTwoHists
from lsst.analysis.tools.vectorActions.selectors import VectorSelector
from lsst.analysis.tools.vectorActions.vectorActions import LoadVector, MagColumnNanoJansky

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
    columnKey = Field(doc="Column key to compute scalars", dtype=str)

    snFluxType = Field(doc="column key for the flux type used in SN selection", dtype=str)

    selector = ConfigurableActionField(doc="Selector to use before computing Scalars", dtype=VectorAction)

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield from super().getInputSchema(**kwargs)

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.sigmaMad = SigmaMadAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.count = CountAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.approxMag = _AproxMedian(  # type: ignore
            vectorKey=self.snFluxType  # type: ignore
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        selection = self.selector(data, **kwargs)  # type: ignore
        if mask is not None:
            mask &= selection
        else:
            mask = selection
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
        self.keyedDataAction = KeyedDataSubsetAction()
        self.keyedDataAction.columnKeys = columns
        self.selectors.flagSelector = CoaddPlotFlagSelector()  # type: ignore
        self.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)  # type: ignore


class ShapeSizeFractionalProcess(KeyedDataAction):
    psfFluxShape = ConfigurableActionField(
        doc="Action to calculate the PSF Shape",
    )
    shapeFracDif = ConfigurableActionField(
        doc="Action which adds shapes to the KeyedData",
    )
    objectSelector = ConfigurableActionField(doc="Action to select which objects should be considered")
    highSNRSelector = ConfigurableActionField(
        doc="Selector action add high SNR stars vector to output",
    )
    lowSNRSelector = ConfigurableActionField(
        doc="Selector action add low SNR stars vector to output",
    )
    calculatorActions = ConfigurableActionStructField(
        doc="Actions which generate KeyedData of scalars",
    )
    magAction = ConfigurableActionField(
        doc="Action to generate the magnitude column"
    )

    def setDefaults(self):
        # compute the PSF Shape
        psfShapeName = "{band}_psfShape"
        self.psfFluxShape = AddComputedVector()
        self.psfFluxShape.keyName = psfShapeName
        self.psfFluxShape.action = CalcShapeSize()
        self.psfFluxShape.action.colXx = "{band}_ixxPSF"
        self.psfFluxShape.action.colYy = "{band}_iyyPSF"
        self.psfFluxShape.action.colXy = "{band}_ixyPSF"

        # compute the Difference of shapes
        self.shapeFracDif = AddComputedVector()
        self.shapeFracDif.keyName = "{band}_derivedShape"
        self.shapeFracDif.action = FractionalDifference()
        self.shapeFracDif.action.actionA = CalcShapeSize()
        self.shapeFracDif.action.actionB = LoadVector()
        self.shapeFracDif.action.actionB.vectorKey = psfShapeName

        # Default mag action action
        self.magAction = MagColumnNanoJansky(columnKey="{band}_psfFlux")

        # Setup the selectors
        self.objectSelector = StellarSelector()
        self.highSNRSelector = SnSelector(threshold=2700)
        self.lowSNRSelector = SnSelector(threshold=500)

        self.calculatorActions.highSNStars = ShapeSizeFractionalScalars(  # type: ignore
            columnKey=self.shapeFracDif.keyName, snFluxType=self.highSNRSelector.fluxType
        )
        self.calculatorActions.highSNStars.selector = VectorSelector(vectorKey='starHighSNMask')  # type: ignore
        self.calculatorActions.lowSNStars = ShapeSizeFractionalScalars(  # type: ignore
            columnKey=self.shapeFracDif.keyName, snFluxType=self.lowSNRSelector.fluxType
        )
        self.calculatorActions.lowSNStars.selector = VectorSelector(vectorKey='starLowSNMask')  # type: ignore

        for action in self.calculatorActions:  # type: ignore
            action.setDefaults()

        super().setDefaults()

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        yield from self.psfFluxShape.getInputSchema(**kwargs)  # type: ignore
        yield from self.objectSelector.getInputSchema(**kwargs)  # type: ignore
        yield from self.highSNRSelector.getInputSchema(**kwargs)  # type: ignore
        yield from self.lowSNRSelector.getInputSchema(**kwargs)  # type: ignore
        yield from self.shapeFracDif.getInputSchema(**kwargs)  # type: ignore
        for action in self.calculatorActions:  # type: ignore
            yield from action.getInputSchema(**kwargs)

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        data = self.psfFluxShape(data, **kwargs)  # type: ignore
        data = self.shapeFracDif(data, **kwargs)  # type: ignore

        objectMask = self.objectSelector(data, **kwargs)  # type: ignore
        highSNRMask = self.highSNRSelector(data, **(kwargs | {"mask": objectMask}))  # type: ignore
        lowSNRMask = self.lowSNRSelector(data, **kwargs | {"mask": objectMask})  # type: ignore

        data['starHighSNMask'] = highSNRMask
        data['starLowSNMask'] = lowSNRMask

        mags = self.magAction(data, **kwargs)  # type: ignore

        for name, action in self.calculatorActions.items():  # type: ignore
            for key, value in action(data, **kwargs).items():
                prefix = f"{band}_" if (band := kwargs.get("band")) else ""
                newKey = f"{prefix}{name}_{key}"
                results[newKey] = value

        results["starHighSNMask"] = highSNRMask[objectMask]
        results["starLowSNMask"] = lowSNRMask[objectMask]
        results['xsStars'] = mags[objectMask]
        shapeDiff = cast(Vector, data[cast(str, self.shapeFracDif.keyName.format(**kwargs))])  # type: ignore
        results["ysStars"] = shapeDiff[objectMask]
        return results


class ShapeSizeFractionalPostProcessPlot(ScatterPlotWithTwoHists):
    def setDefaults(self):
        self.plotTypes = ["stars"]
        self.xAxisLabel = "PSF Magnitude (mag)"
        self.yAxisLabel = "Fractional size residuals (S/S_PSF - 1)"
        self.magLabel = "PSF Magnitude (mag)"
        self.highThreshold = 2700
        self.lowThreshold = 500


class ShapeSizeFractionalPostProcessMetric(MetricAction):
    units = DictField(
        doc="Mapping of scalar key to astropy unit string",
        keytype=str,
        itemtype=str,
        default={
            "{band}_highSNStars_median": "pixel",
            "{band}_highSNStars_sigmaMad": "pixel",
            "{band}_highSNStars_count": "count",
            "{band}_lowSNStars_median": "pixel",
            "{band}_lowSNStars_sigmaMad": "pixel",
            "{band}_lowSNStars_count": "count",
        },
    )
    newNames = DictField(
        doc="Mapping of key to new name if needed prior to creating metric",
        keytype=str,
        itemtype=str,
        default={},
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
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
