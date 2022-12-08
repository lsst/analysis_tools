from lsst.analysis.tools.actions.scalar import CountUniqueAction
from lsst.analysis.tools.actions.vector import DownselectVector, LoadVector, ThresholdSelector, VectorSelector
from lsst.analysis.tools.interfaces import AnalysisTool

__all__ = ("BaseSources",)


class BaseSources(AnalysisTool):
    """Base class for counting associated and unassociated sources."""

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.loadVectorSources = LoadVector()
        self.process.buildActions.associatedVectorSelector = ThresholdSelector()
        self.process.buildActions.unassociatedVectorSelector = ThresholdSelector()

        # assign keys for PSF and AP Flux
        self.process.buildActions.loadVectorSources.vectorKey = "nDiaSources"
        self.process.buildActions.associatedVectorSelector.vectorKey = "nDiaSources"
        self.process.buildActions.associatedVectorSelector.op = "gt"
        self.process.buildActions.associatedVectorSelector.threshold = 1.0
        self.process.buildActions.unassociatedVectorSelector.vectorKey = "nDiaSources"
        self.process.buildActions.unassociatedVectorSelector.op = "le"
        self.process.buildActions.unassociatedVectorSelector.threshold = 1.0

        self.process.filterActions.allSources = VectorSelector(vectorKey="nDiaSources")
        self.process.filterActions.associatedVector = DownselectVector(
            vectorKey="nDiaSources", selector=self.process.buildActions.associatedVectorSelector
        )
        self.process.filterActions.unassociatedVector = DownselectVector(
            vectorKey="nDiaSources", selector=self.process.buildActions.unassociatedVectorSelector
        )

        self.process.buildActions.uniqueSources = CountUniqueAction(vectorKey="nDiaSources")
