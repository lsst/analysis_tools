# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregate metrics from the standard configuration of
# lsst.analysis.tools.tasks.AssociatedSourcesTractAnalysisTask.

from lsst.analysis.tools.atools import WholeSkyPlotTool

# Keys with band:
keysWithBand = [
    "coadd_depth_summary_metrics_depth_above_threshold_1_{band}_mean",
    "coadd_depth_summary_metrics_depth_above_threshold_3_{band}_mean",
    "coadd_depth_summary_metrics_depth_above_threshold_5_{band}_mean",
    "coadd_depth_summary_metrics_depth_above_threshold_12_{band}_mean",
]

for key in keysWithBand:
    atoolName = key.replace("_{band}", "")
    setattr(config.atools, atoolName, WholeSkyPlotTool)

    atool = getattr(config.atools, atoolName)
    setattr(atool, "metric", key)

    plot = getattr(getattr(atool, "produce"), "plot")
    setattr(plot, "showOutliers", False)
    setattr(plot, "showNaNs", False)
    setattr(plot, "labelTracts", True)
    setattr(plot, "colorBarMin", 0.0)
    setattr(plot, "colorBarMax", 100.0)
    setattr(plot, "colorBarRange", 1.0)

config.addOutputNamePrefix = True
