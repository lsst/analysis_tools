
__all__ = (
    "DemoMetric",
)


from ..interfaces import AnalysisMetric
from lsst.analysis.tools.actions.scalar import MedianAction
from lsst.analysis.tools.actions.vector import SnSelector


class DemoMetric(AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()

        # select on high signal to noise obejcts
        # add in a signal to noise selector
        self.prep.selectors.snSelector = SnSelector()

        # set what key the selector should use when deciding SNR
        self.prep.selectors.snSelector.fluxType = "psFlux"

        # select what threshold value is desireable for the selector
        self.prep.selectors.snSelector.threshold = 10

        # the final name in the qualification is used as a key to insert
        # the calculation into KeyedData
        self.process.calculateActions.medianValueName = MedianAction(vectorKey="psFlux")

        # tell the metic what the units are for the quantity
        self.produce.units = {"medianValueName": "Jy"}

        # Rename the quanity prior to producing the Metric
        # (useful for resuable workflows that set a name toward
        # the end of computation)
        self.produce.newNames = {"medianValueName": "DemoMetric"}
