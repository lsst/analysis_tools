# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregates metrics from objectTableCore analysis.

from lsst.analysis.tools.atools import WholeSkyPlotTool

# Keys with `band`
keysWithBand = [
    "e1Diff_{band}_highSNStars_median",
    "e1Diff_{band}_highSNStars_sigmaMad",
    "e1Diff_{band}_lowSNStars_median",
    "e1Diff_{band}_lowSNStars_sigmaMad",
    "e2Diff_{band}_highSNStars_median",
    "e2Diff_{band}_highSNStars_sigmaMad",
    "e2Diff_{band}_lowSNStars_median",
    "e2Diff_{band}_lowSNStars_sigmaMad",
    "psfCModelScatter_{band}_psf_cModel_diff_median",
    "psfCModelScatter_{band}_psf_cModel_diff_sigmaMad",
    "psfCModelScatter_{band}_psf_cModel_diff_mean",
    "shapeSizeFractionalDiff_{band}_highSNStars_median",
    "shapeSizeFractionalDiff_{band}_highSNStars_sigmaMad",
    "shapeSizeFractionalDiff_{band}_lowSNStars_median",
    "shapeSizeFractionalDiff_{band}_lowSNStars_sigmaMad",
    "skyFluxStatisticMetric_{band}_medianSky",
    "skyFluxStatisticMetric_{band}_meanSky",
    "skyFluxStatisticMetric_{band}_stdevSky",
    "skyFluxStatisticMetric_{band}_sigmaMADSky",
]
if hasattr(parameters, "objectTableCoreWholeSkyPlotKeysWithBand"):
    keysWithBand = parameters.objectTableCoreWholeSkyPlotKeysWithBand

for keyWithBand in keysWithBand:
    atoolName = keyWithBand.replace("_{band}", "")
    setattr(config.atools, atoolName, WholeSkyPlotTool)
    atool = getattr(config.atools, atoolName)
    setattr(atool, "metric", keyWithBand)

# Keys without `band`
keys = [
    "wPerpPSF_wPerp_psfFlux_median",
    "wPerpPSF_wPerp_psfFlux_sigmaMAD",
    "yPerpPSF_yPerp_psfFlux_median",
    "yPerpPSF_yPerp_psfFlux_sigmaMAD",
]
if hasattr(parameters, "objectTableCoreWholeSkyPlotKeys"):
    keys = parameters.objectTableCoreWholeSkyPlotKeys

for atoolName in keys:
    setattr(config.atools, atoolName, WholeSkyPlotTool)
    atool = getattr(config.atools, atoolName)
    setattr(atool, "metric", atoolName)
    setattr(atool, "parameterizedBand", False)

config.addOutputNamePrefix = True
