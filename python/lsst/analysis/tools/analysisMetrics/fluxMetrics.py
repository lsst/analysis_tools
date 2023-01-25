__all__ = ("CentralTendency",)

from lsst.analysis.tools.actions.scalar import MeanAction, MedianAction
from lsst.analysis.tools.interfaces import AnalysisMetric


class CentralTendency(AnalysisMetric):
    """Metric for measuring mean and median of psf, ap,
    and total flux.
    """

    def setDefaults(self):
        super().setDefaults()

        self.process.calculateActions.psFluxMedianMetric = MedianAction(vectorKey="psFlux")
        self.process.calculateActions.apFluxMedianMetric = MedianAction(vectorKey="apFlux")
        self.process.calculateActions.totFluxMedianMetric = MedianAction(vectorKey="totFlux")

        self.process.calculateActions.psFluxMeanMetric = MeanAction(vectorKey="psFlux")
        self.process.calculateActions.apFluxMeanMetric = MeanAction(vectorKey="apFlux")
        self.process.calculateActions.totFluxMeanMetric = MeanAction(vectorKey="totFlux")

        self.produce.units = {
            "psFluxMeanMetric": "flx",
            "apFluxMeanMetric": "flx",
            "totFluxMeanMetric": "flx",
            "psFluxMedianMetric": "flx",
            "apFluxMedianMetric": "flx",
            "totFluxMedianMetric": "flx",
        }
