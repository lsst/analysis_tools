# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregates metrics from objectParentTableCore analysis.

from lsst.analysis.tools.atools import WholeSkyPlotTool

# Keys without `band`
keys = [
    "skippedDeblenderMetrics_numSkippedPeaks",
    "skippedDeblenderMetrics_numSkippedBlends",
    "skippedDeblenderMetrics_numBlendParentTooBig",
    "skippedDeblenderMetrics_numBlendTooManyPeaks",
    "skippedDeblenderMetrics_numBlendTooManyMasked",
]
if hasattr(parameters, "objectParentTableCoreWholeSkyPlotKeys"):
    keys = parameters.objectParentTableCoreWholeSkyPlotKeys

for atoolName in keys:
    setattr(config.atools, atoolName, WholeSkyPlotTool)
    atool = getattr(config.atools, atoolName)
    setattr(atool, "metric", atoolName)
    setattr(atool, "parameterizedBand", False)

config.addOutputNamePrefix = True
