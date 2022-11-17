from lsst.analysis.tools.interfaces import AnalysisTool
from lsst.analysis.tools.actions.vector import LoadVector, DivideVector


class DemoPsfApRatioBaseClass(AnalysisTool):
    """Base class for scatter plots of PSF residuals.
    This is shared by size and ellipticity plots.
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.loadVectorPsf = LoadVector()
        self.process.filterActions.loadVectorAp = LoadVector()

        self.process.filterActions.loadVectorPsf.vectorKey = "psFlux"
        self.process.filterActions.loadVectorAp.vectorKey = "apFlux"

        # the final name in the qualification is used as a key to insert
        # the calculation into KeyedData
        self.process.calculateActions.fluxRatioMetric = DivideVector(
            actionA=self.process.filterActions.loadVectorPsf,
            actionB=self.process.filterActions.loadVectorAp)
