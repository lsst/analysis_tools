# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregates metrics from objectTableCoreRefCatMatch analysis.

from lsst.analysis.tools.atools import WholeSkyPlotTool

# Keys with `band`
keysWithBand = [
    "astromDiffMetrics_{band}_AA1_RA_coadd",
    "astromDiffMetrics_{band}_AA1_sigmaMad_RA_coadd",
    "astromDiffMetrics_{band}_AA1_Dec_coadd",
    "astromDiffMetrics_{band}_AA1_sigmaMad_Dec_coadd",
    "astromDiffMetrics_{band}_AA1_tot_coadd",
    "astromDiffMetrics_{band}_AA1_sigmaMad_tot_coadd",
]
if hasattr(parameters, "objectTableCoreRefCatMatchWholeSkyPlotKeysWithBand"):
    keysWithBand = parameters.objectTableCoreRefCatMatchWholeSkyPlotKeysWithBand

for keyWithBand in keysWithBand:
    atoolName = keyWithBand.replace("_{band}", "")
    setattr(config.atools, atoolName, WholeSkyPlotTool)
    atool = getattr(config.atools, atoolName)
    setattr(atool, "metric", keyWithBand)

config.addOutputNamePrefix = True
