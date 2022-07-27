from __future__ import annotations

from typing import cast

import numpy as np
from lsst.analysis.tools.actions.vector import LoadVector, MagColumnNanoJansky, VectorSelector
from lsst.pex.config import Field
from lsst.pipe.tasks.configurableActions import ConfigurableActionField, ConfigurableActionStructField

from ..actions.keyedData import AddComputedVector, KeyedScalars
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction
from ..actions.scalar import CountAction, MedianAction, SigmaMadAction
from ..actions.vector import (
    CalcE,
    CalcEDiff,
    CalcShapeSize,
    CoaddPlotFlagSelector,
    DownselectVector,
    FractionalDifference,
    SnSelector,
    StarSelector,
)
from ..interfaces import (
    AnalysisAction,
    AnalysisTool,
    KeyedData,
    KeyedDataAction,
    KeyedDataSchema,
    Scalar,
    ScalarAction,
    Vector,
    VectorAction,
)


class _ApproxMedian(ScalarAction):
    vectorKey = Field[str](doc="Key for the vector to perform action on", optional=False)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        mask = self.getMask(**kwargs)
        value = np.sort(data[self.vectorKey.format(**kwargs)][mask])  # type: ignore
        x = int(len(value) / 10)
        return np.nanmedian(value[-x:])


class ShapeSizeFractionalScalars(KeyedScalars):
    vectorKey = Field[str](doc="Column key to compute scalars")

    snFluxType = Field[str](doc="column key for the flux type used in SN selection")

    selector = ConfigurableActionField[VectorAction](doc="Selector to use before computing Scalars")

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.sigmaMad = SigmaMadAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.count = CountAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.approxMag = _ApproxMedian(vectorKey=self.snFluxType)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        selection = self.selector(data, **kwargs)
        if mask is not None:
            mask &= selection
        else:
            mask = selection
        return super().__call__(data, **kwargs | dict(mask=mask))


class ShapeSizeFractionalProcessLegacy(KeyedDataAction):
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
        self.magAction = MagColumnNanoJansky(vectorKey="{band}_psfFlux")

        # Setup the selectors
        self.objectSelector = StarSelector()
        self.highSNRSelector = SnSelector(threshold=2700)
        self.lowSNRSelector = SnSelector(threshold=500)

        self.calculatorActions.highSNStars = ShapeSizeFractionalScalars(
            vectorKey=self.shapeFracDif.keyName, snFluxType=self.highSNRSelector.fluxType
        )
        highSNSelector = VectorSelector(vectorKey="starHighSNMask")
        self.calculatorActions.highSNStars.selector = highSNSelector
        self.calculatorActions.lowSNStars = ShapeSizeFractionalScalars(
            vectorKey=self.shapeFracDif.keyName, snFluxType=self.lowSNRSelector.fluxType
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
                outputs = action.getOutputSchema()
                if outputs is not None:
                    yield from outputs
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


class BasePsfResidualMixin(AnalysisTool):
    """Shared configuration for `prep` and `process` stages of PSF residuals.

    This is a mixin class used by `BasePsfResidualScatterPlot` and
    `BasePsfResidualMetric` to share common default configuration.
    """

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.snSelector = SnSelector(fluxType="{band}_psfFlux", threshold=100)

        self.process.buildActions.mags = MagColumnNanoJansky(vectorKey="{band}_psfFlux")
        self.process.buildActions.fracDiff = FractionalDifference(
            actionA=CalcShapeSize(colXx="{band}_ixx", colYy="{band}_iyy", colXy="{band}_ixy"),
            actionB=CalcShapeSize(colXx="{band}_ixxPSF", colYy="{band}_iyyPSF", colXy="{band}_ixyPSF"),
        )
        # Define an eDiff action and let e1Diff and e2Diff differ only in
        # component.
        self.process.buildActions.eDiff = CalcEDiff(
            colA=CalcE(colXx="{band}_ixx", colYy="{band}_iyy", colXy="{band}_ixy"),
            colB=CalcE(colXx="{band}_ixxPSF", colYy="{band}_iyyPSF", colXy="{band}_ixyPSF"),
        )
        self.process.buildActions.e1Diff = self.process.buildActions.eDiff
        self.process.buildActions.e1Diff.component = "1"
        self.process.buildActions.e2Diff = self.process.buildActions.eDiff
        self.process.buildActions.e2Diff.component = "2"
        # pre-compute a stellar selector mask so it can be used in the filter
        # actions while only being computed once, alternatively the stellar
        # selector could be calculated and applied twice in the filter stage
        self.process.buildActions.starSelector = StarSelector()

        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags", selector=VectorSelector(vectorKey="starSelector")
        )
        # downselect the psfFlux as well
        self.process.filterActions.psfFlux = DownselectVector(
            vectorKey="{band}_psfFlux", selector=VectorSelector(vectorKey="starSelector")
        )
        self.process.filterActions.psfFluxErr = DownselectVector(
            vectorKey="{band}_psfFluxErr", selector=VectorSelector(vectorKey="starSelector")
        )

        self.process.calculateActions.stars = ScatterPlotStatsAction(
            vectorKey="yStars",
        )
        # use the downselected psfFlux
        self.process.calculateActions.stars.highSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.lowSNSelector.fluxType = "psfFlux"
        self.process.calculateActions.stars.fluxType = "psfFlux"
