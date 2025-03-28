# Configuration for standard instances of
# lsst.analysis.tools.tasks.AssociatedSourcesTractAnalysisTask

from lsst.analysis.tools.atools import *

config.atools.stellarAstrometricSelfRepeatabilityRA = AstrometricRepeatability
config.atools.stellarAstrometricSelfRepeatabilityRA.coordinate = "RA"
config.atools.stellarAstrometricSelfRepeatabilityRA.level = 2
config.atools.stellarAstrometricSelfRepeatabilityDec = AstrometricRepeatability
config.atools.stellarAstrometricSelfRepeatabilityDec.coordinate = "Dec"
config.atools.stellarAstrometricSelfRepeatabilityDec.level = 2
config.atools.stellarPhotometricRepeatability = StellarPhotometricRepeatability
config.atools.stellarPhotometricResiduals = StellarPhotometricResidualsFocalPlane
config.atools.stellarPhotometricRepeatabilityCalib = StellarPhotometricRepeatability
config.atools.stellarPhotometricRepeatabilityCalib.fluxType = "calibFlux"
config.atools.stellarPhotometricResidualsCalib = StellarPhotometricResidualsFocalPlane
config.atools.stellarPhotometricResidualsCalib.fluxType = "calibFlux"
config.atools.stellarAstrometricResidualsRA = StellarAstrometricResidualsRAFocalPlanePlot
config.atools.stellarAstrometricResidualsDec = StellarAstrometricResidualsDecFocalPlanePlot
config.atools.stellarAstrometricResidualStdDevRA = StellarAstrometricResidualStdDevRAFocalPlanePlot
config.atools.stellarAstrometricResidualStdDevDec = StellarAstrometricResidualStdDevDecFocalPlanePlot
config.atools.stellarAstrometricRepeatability1 = AstrometricRelativeRepeatability
config.atools.stellarAstrometricRepeatability1.xValue = 1
config.atools.stellarAstrometricRepeatability1.process.calculateActions.rms.annulus = 5
config.atools.stellarAstrometricRepeatability2 = AstrometricRelativeRepeatability
config.atools.stellarAstrometricRepeatability2.xValue = 2
config.atools.stellarAstrometricRepeatability2.process.calculateActions.rms.annulus = 20
config.atools.stellarAstrometricRepeatability3 = AstrometricRelativeRepeatability
config.atools.stellarAstrometricRepeatability3.xValue = 3
config.atools.stellarAstrometricRepeatability3.process.calculateActions.rms.annulus = 200
config.atools.stellarAstrometricRepeatability3.process.calculateActions.rms.threshAD = 30
config.addOutputNamePrefix = True
