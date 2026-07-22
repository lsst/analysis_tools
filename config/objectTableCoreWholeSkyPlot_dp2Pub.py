# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregates metrics from objectTableCore analysis.

from lsst.analysis.tools.atools import WholeSkyPlotTool

# Keys with `band`
keysWithBand = [
    "psfCModelScatter_{band}_psf_cModel_diff_median",
    "shapeSizeFractionalDiff_{band}_highSNStars_median",
    "psfSersicScatter_{band}_psf_sersic_diff_median",
    "fiveSigmaDepth_{band}_median5sigmaDepth",
]
if hasattr(parameters, "objectTableCoreWholeSkyPlotKeysWithBand"):
    keysWithBand = parameters.objectTableCoreWholeSkyPlotKeysWithBand

for keyWithBand in keysWithBand:
    atoolName = keyWithBand.replace("_{band}", "")
    setattr(config.atools, atoolName, WholeSkyPlotTool)
    atool = getattr(config.atools, atoolName)
    setattr(atool, "metric", keyWithBand)
    setattr(atool, "publicationStyle", True)
    if "Depth" not in atoolName:
        setattr(atool, "fixAroundZero", True)

config.addOutputNamePrefix = True
