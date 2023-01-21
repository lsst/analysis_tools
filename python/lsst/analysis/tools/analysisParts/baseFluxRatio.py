__all__ = ("BasePsfApRatio",)

from lsst.analysis.tools.actions.vector import DivideVector, LoadVector
from lsst.analysis.tools.interfaces import AnalysisTool


class BasePsfApRatio(AnalysisTool):
    """Base class for plots or metrics which use PSF/Aperture Ratios."""

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.loadVectorPsf = LoadVector()
        self.process.buildActions.loadVectorAp = LoadVector()

        # assign keys for PSF and AP Flux
        self.process.buildActions.loadVectorPsf.vectorKey = "psFlux"
        self.process.buildActions.loadVectorAp.vectorKey = "apFlux"

        self.process.calculateActions.fluxRatio = DivideVector(
            actionA=self.process.buildActions.loadVectorPsf, actionB=self.process.buildActions.loadVectorAp
        )
