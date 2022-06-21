from __future__ import annotations

import numpy as np

from typing import cast
from lsst.analysis.tools.vectorActions.selectors import VectorSelector
from lsst.analysis.tools.vectorActions.vectorActions import LoadVector, MagColumnNanoJansky

from lsst.pex.config import Field
from lsst.pipe.tasks.configurableActions import ConfigurableActionField, ConfigurableActionStructField

from ..interfaces import (
    AnalysisAction,
    KeyedData,
    KeyedDataAction,
    Scalar,
    ScalarAction,
    VectorAction,
    KeyedDataSchema,
    Vector,
)
from ..scalarActions import CountAction, MedianAction, SigmaMadAction
from ..keyedDataActions import AddComputedVector, KeyedScalars
from ..vectorActions import SnSelector, StellarSelector, FractionalDifference
from ..vectorActions.calcShapeSize import CalcShapeSize


class _AproxMedian(ScalarAction):
    vectorKey = Field[str](doc="Key for the vector to perform action on", optional=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(data[self.vectorKey.format(**kwargs)][mask])  # type: ignore
        x = int(len(value) / 10)
        return np.nanmedian(value[-x:])


class ShapeSizeFractionalScalars(KeyedScalars):
    columnKey = Field[str](doc="Column key to compute scalars")

    snFluxType = Field[str](doc="column key for the flux type used in SN selection")

    selector = ConfigurableActionField[VectorAction](doc="Selector to use before computing Scalars")

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.sigmaMad = SigmaMadAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.count = CountAction(colKey=self.columnKey)  # type: ignore
        self.scalarActions.approxMag = _AproxMedian(  # type: ignore
            vectorKey=self.snFluxType
        )

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        selection = self.selector(data, **kwargs)
        if mask is not None:
            mask &= selection
        else:
            mask = selection
        return super().__call__(data, **kwargs | dict(mask=mask))


class ShapeSizeFractionalProcess(KeyedDataAction):
    psfFluxShape = ConfigurableActionField[AnalysisAction](
        doc="Action to calculate the PSF Shape",
    )
    shapeFracDif = ConfigurableActionField[AnalysisAction](
        doc="Action which adds shapes to the KeyedData",
    )
    objectSelector = ConfigurableActionField[AnalysisAction](
        doc="Action to select which objects should be considered"
    )
    highSNRSelector = ConfigurableActionField[AnalysisAction](
        doc="Selector action add high SNR stars vector to output",
    )
    lowSNRSelector = ConfigurableActionField[AnalysisAction](
        doc="Selector action add low SNR stars vector to output",
    )
    calculatorActions = ConfigurableActionStructField[KeyedDataAction](
        doc="Actions which generate KeyedData of scalars",
    )
    magAction = ConfigurableActionField(doc="Action to generate the magnitude column")

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

        self.calculatorActions.highSNStars = ShapeSizeFractionalScalars(
            columnKey=self.shapeFracDif.keyName, snFluxType=self.highSNRSelector.fluxType
        )
        highSNSelector = VectorSelector(vectorKey="starHighSNMask")
        self.calculatorActions.highSNStars.selector = highSNSelector
        self.calculatorActions.lowSNStars = ShapeSizeFractionalScalars(
            columnKey=self.shapeFracDif.keyName, snFluxType=self.lowSNRSelector.fluxType
        )
        self.calculatorActions.lowSNStars.selector = VectorSelector(vectorKey="starLowSNMask")

        for action in self.calculatorActions:
            action.setDefaults()

        super().setDefaults()

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.psfFluxShape.getInputSchema()
        yield from self.objectSelector.getInputSchema()
        yield from self.highSNRSelector.getInputSchema()
        yield from self.lowSNRSelector.getInputSchema()
        yield from self.shapeFracDif.getInputSchema()
        for action in self.calculatorActions:
            yield from action.getInputSchema()

    def getOutputSchema(self) -> KeyedDataSchema:
        results = (
            ("starHighSNMask", Vector),
            ("starLowSNMask", Vector),
            ("xsStars", Vector),
            ("ysStars", Vector),
        )
        for action in self.calculatorActions:
            if isinstance(action, KeyedDataAction):
                yield from action.getOutputSchema()
        yield from results

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        results = {}
        data = self.psfFluxShape(data, **kwargs)
        data = self.shapeFracDif(data, **kwargs)

        objectMask = self.objectSelector(data, **kwargs)
        highSNRMask = self.highSNRSelector(data, **(kwargs | {"mask": objectMask}))
        lowSNRMask = self.lowSNRSelector(data, **kwargs | {"mask": objectMask})

        data["starHighSNMask"] = highSNRMask
        data["starLowSNMask"] = lowSNRMask

        mags = self.magAction(data, **kwargs)  # type: ignore

        for name, action in self.calculatorActions.items():  # type: ignore
            for key, value in action(data, **kwargs).items():
                prefix = f"{band}_" if (band := kwargs.get("band")) else ""
                newKey = f"{prefix}{name}_{key}"
                results[newKey] = value

        results["starHighSNMask"] = highSNRMask[objectMask]
        results["starLowSNMask"] = lowSNRMask[objectMask]
        results["xsStars"] = mags[objectMask]
        shapeDiff = cast(Vector, data[cast(str, self.shapeFracDif.keyName.format(**kwargs))])  # type: ignore
        results["ysStars"] = shapeDiff[objectMask]
        return results
