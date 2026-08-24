# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregates metrics from objectTableCoreRefCatMatch analysis.

from lsst.analysis.tools.atools import WholeSkyPlotTool

# Keys with `band`
keysWithBand = [
    "astromDiffRAScatterPlot_{band}_ref_ra_offset_coadd",
    "astromDiffRAScatterPlot_{band}_ref_ra_offset_sigmaMad_coadd",
    "astromDiffDecScatterPlot_{band}_ref_dec_offset_coadd",
    "astromDiffDecScatterPlot_{band}_ref_dec_offset_sigmaMad_coadd",
]
if hasattr(parameters, "objectTableCoreRefCatMatchWholeSkyPlotKeysWithBand"):
    keysWithBand = parameters.objectTableCoreRefCatMatchWholeSkyPlotKeysWithBand

for keyWithBand in keysWithBand:
    atoolName = keyWithBand.replace("_{band}", "")
    setattr(config.atools, atoolName, WholeSkyPlotTool)
    atool = getattr(config.atools, atoolName)
    setattr(atool, "metric", keyWithBand)
    setattr(atool, "publicationStyle", True)
    if "sigmaMad" not in atoolName:
        setattr(atool, "fixAroundZero", True)

config.addOutputNamePrefix = True
