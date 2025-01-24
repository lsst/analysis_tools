# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregate metrics from the standard configuration of
# lsst.analysis.tools.tasks.AssociatedSourcesTractAnalysisTask.

from lsst.analysis.tools.atools import *

config.atools.wholeSkyMetric = WholeSkyPlotTool
config.atools.wholeSkyMetric.plotKeys = []
config.atools.wholeSkyMetric.keysWithBand = [
    "stellarAstrometricRepeatability1_{band}_AM1",
    "stellarAstrometricRepeatability1_{band}_AF1",
    "stellarAstrometricRepeatability1_{band}_AD1",
    "stellarAstrometricRepeatability2_{band}_AM2",
    "stellarAstrometricRepeatability2_{band}_AF2",
    "stellarAstrometricRepeatability2_{band}_AD2",
    "stellarAstrometricRepeatability3_{band}_AM3",
    "stellarAstrometricRepeatability3_{band}_AF3",
    "stellarAstrometricRepeatability3_{band}_AD3",
    "stellarAstrometricSelfRepeatabilityDec_{band}_dmL2AstroErr_Dec",
    "stellarAstrometricSelfRepeatabilityRA_{band}_dmL2AstroErr_RA",
    "stellarPhotometricRepeatability_{band}_stellarPhotRepeatStdev",
    "stellarPhotometricRepeatability_{band}_stellarPhotRepeatOutlierFraction",
    "stellarPhotometricResiduals_{band}_photResidTractSigmaMad",
    "stellarPhotometricResiduals_{band}_photResidTractStdev",
    "stellarPhotometricResiduals_{band}_photResidTractMedian",
]
