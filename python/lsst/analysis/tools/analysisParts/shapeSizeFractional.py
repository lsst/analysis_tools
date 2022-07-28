from __future__ import annotations

from lsst.pex.config import Field
from lsst.pipe.tasks.configurableActions import ConfigurableActionField

from ..actions.keyedData import KeyedScalars
from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction
from ..actions.scalar import CountAction, MedianAction, SigmaMadAction
from ..actions.vector import (
    CalcE,
    CalcEDiff,
    CalcShapeSize,
    CoaddPlotFlagSelector,
    DownselectVector,
    FractionalDifference,
    MagColumnNanoJansky,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from ..interfaces import AnalysisTool, KeyedData, VectorAction


class ShapeSizeFractionalScalars(KeyedScalars):
    vectorKey = Field[str](doc="Column key to compute scalars")

    snFluxType = Field[str](doc="column key for the flux type used in SN selection")

    selector = ConfigurableActionField[VectorAction](doc="Selector to use before computing Scalars")

    def setDefaults(self):
        super().setDefaults()
        self.scalarActions.median = MedianAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.sigmaMad = SigmaMadAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.count = CountAction(vectorKey=self.vectorKey)  # type: ignore
        self.scalarActions.approxMag = MedianAction(vectorKey=self.snFluxType)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        mask = kwargs.get("mask")
        selection = self.selector(data, **kwargs)
        if mask is not None:
            mask &= selection
        else:
            mask = selection
        return super().__call__(data, **kwargs | dict(mask=mask))


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
