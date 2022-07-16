from __future__ import annotations

__all__ = ("Ap12_PSF_skyPlot",)

from ..vectorActions.selectors import VisitPlotFlagSelector, SnSelector, StellarSelector
from ..plotActions.scatterplotWithTwoHists import ScatterPlotWithTwoHists, ScatterPlotStatsAction
from ..plotActions.colorColorFitPlot import ColorColorFitPlot
from ..plotActions.skyPlot import SkyPlot
from ..vectorActions.vectorActions import (
    MagDiff,
    MagColumnNanoJansky,
    DownselectVector,
    VectorSelector,
    FractionalDifference,
    LoadVector,
)
from ..vectorActions.calcShapeSize import CalcShapeSize

from ..scalarActions.scalarActions import ApproxFloor
from ..keyedDataActions.stellarLocusFit import StellarLocusFitAction

from ..interfaces import AnalysisPlot


class Ap12_PSF_skyPlot(AnalysisPlot):
    band: str = "{band}"

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "psfFlux"
        self.prep.selectors.snSelector.threshold = 100

        self.prep.selectors.starSelector = StellarSelector()
        self.prep.selectors.starSelector.columnKey = "extendedness"

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "psfFlux"

        self.process.buildActions.zStars = MagDiff()
        self.process.buildActions.zStars.col1 = f"ap12Flux"
        self.process.buildActions.zStars.col2 = f"psfFlux"

        self.post_process = SkyPlot()
        self.post_process.plotTypes = ["stars"]
        self.post_process.plotName = f"ap12-psf_{self.band}"
        self.post_process.xAxisLabel = "R.A."
        self.post_process.yAxisLabel = "Declination"
        self.post_process.zAxisLabel = "Ap 12 - PSF [mag]"
        self.post_process.plotOutlines = False
