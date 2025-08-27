# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregate metrics from the standard configuration of
# lsst.analysis.tools.tasks.AssociatedSourcesTractAnalysisTask.

from lsst.analysis.tools.atools import WholeSkyPlotTool

config.atools.stellarAstrometricRepeatability1_AM1 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability1_AM1.metric = "stellarAstrometricRepeatability1_{band}_AM1"
config.atools.stellarAstrometricRepeatability1_AF1 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability1_AF1.metric = "stellarAstrometricRepeatability1_{band}_AF1"
config.atools.stellarAstrometricRepeatability1_AD1 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability1_AD1.metric = "stellarAstrometricRepeatability1_{band}_AD1"
config.atools.stellarAstrometricRepeatability2_AM2 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability2_AM2.metric = "stellarAstrometricRepeatability2_{band}_AM2"
config.atools.stellarAstrometricRepeatability2_AF2 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability2_AF2.metric = "stellarAstrometricRepeatability2_{band}_AF2"
config.atools.stellarAstrometricRepeatability2_AD2 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability2_AD2.metric = "stellarAstrometricRepeatability2_{band}_AD2"
config.atools.stellarAstrometricRepeatability3_AM3 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability3_AM3.metric = "stellarAstrometricRepeatability3_{band}_AM3"
config.atools.stellarAstrometricRepeatability3_AF3 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability3_AF3.metric = "stellarAstrometricRepeatability3_{band}_AF3"
config.atools.stellarAstrometricRepeatability3_AD3 = WholeSkyPlotTool
config.atools.stellarAstrometricRepeatability3_AD3.metric = "stellarAstrometricRepeatability3_{band}_AD3"
config.atools.stellarAstrometricSelfRepeatabilityDec_dmL2AstroErr_Dec = WholeSkyPlotTool
config.atools.stellarAstrometricSelfRepeatabilityDec_dmL2AstroErr_Dec.metric = (
    "stellarAstrometricSelfRepeatabilityDec_{band}_dmL2AstroErr_Dec"
)
config.atools.stellarAstrometricSelfRepeatabilityRA_dmL2AstroErr_RA = WholeSkyPlotTool
config.atools.stellarAstrometricSelfRepeatabilityRA_dmL2AstroErr_RA.metric = (
    "stellarAstrometricSelfRepeatabilityRA_{band}_dmL2AstroErr_RA"
)
config.atools.stellarPhotometricRepeatability_stellarPhotRepeatStdev = WholeSkyPlotTool
config.atools.stellarPhotometricRepeatability_stellarPhotRepeatStdev.metric = (
    "stellarPhotometricRepeatability_{band}_stellarPhotRepeatStdev"
)
config.atools.stellarPhotometricRepeatability_stellarPhotRepeatOutlierFraction = WholeSkyPlotTool
config.atools.stellarPhotometricRepeatability_stellarPhotRepeatOutlierFraction.metric = (
    "stellarPhotometricRepeatability_{band}_stellarPhotRepeatOutlierFraction"
)
config.atools.stellarPhotometricResiduals_photResidTractSigmaMad = WholeSkyPlotTool
config.atools.stellarPhotometricResiduals_photResidTractSigmaMad.metric = (
    "stellarPhotometricResiduals_{band}_photResidTractSigmaMad"
)
config.atools.stellarPhotometricResiduals_photResidTractStdev = WholeSkyPlotTool
config.atools.stellarPhotometricResiduals_photResidTractStdev.metric = (
    "stellarPhotometricResiduals_{band}_photResidTractStdev"
)
config.atools.stellarPhotometricResiduals_photResidTractMedian = WholeSkyPlotTool
config.atools.stellarPhotometricResiduals_photResidTractMedian.metric = (
    "stellarPhotometricResiduals_{band}_photResidTractMedian"
)
config.addOutputNamePrefix = True
