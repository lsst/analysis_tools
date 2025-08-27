# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask
# that aggregates metrics from objectTableCore analysis.

from lsst.analysis.tools.atools import WholeSkyPlotTool

# Keys With `band`
config.atools.e1Diff_highSNStars_median = WholeSkyPlotTool
config.atools.e1Diff_highSNStars_median.metric = "e1Diff_{band}_highSNStars_median"
config.atools.e1Diff_highSNStars_median = WholeSkyPlotTool
config.atools.e1Diff_highSNStars_median.metric = "e1Diff_{band}_highSNStars_median"
config.atools.e1Diff_highSNStars_sigmaMad = WholeSkyPlotTool
config.atools.e1Diff_highSNStars_sigmaMad.metric = "e1Diff_{band}_highSNStars_sigmaMad"
config.atools.e1Diff_lowSNStars_median = WholeSkyPlotTool
config.atools.e1Diff_lowSNStars_median.metric = "e1Diff_{band}_lowSNStars_median"
config.atools.e1Diff_lowSNStars_sigmaMad = WholeSkyPlotTool
config.atools.e1Diff_lowSNStars_sigmaMad.metric = "e1Diff_{band}_lowSNStars_sigmaMad"
config.atools.e2Diff_highSNStars_median = WholeSkyPlotTool
config.atools.e2Diff_highSNStars_median.metric = "e2Diff_{band}_highSNStars_median"
config.atools.e2Diff_highSNStars_sigmaMad = WholeSkyPlotTool
config.atools.e2Diff_highSNStars_sigmaMad.metric = "e2Diff_{band}_highSNStars_sigmaMad"
config.atools.e2Diff_lowSNStars_median = WholeSkyPlotTool
config.atools.e2Diff_lowSNStars_median.metric = "e2Diff_{band}_lowSNStars_median"
config.atools.e2Diff_lowSNStars_sigmaMad = WholeSkyPlotTool
config.atools.e2Diff_lowSNStars_sigmaMad.metric = "e2Diff_{band}_lowSNStars_sigmaMad"
config.atools.psfCModelScatter_psf_cModel_diff_median = WholeSkyPlotTool
config.atools.psfCModelScatter_psf_cModel_diff_median.metric = (
    "psfCModelScatter_{band}_psf_cModel_diff_median"
)
config.atools.psfCModelScatter_psf_cModel_diff_sigmaMad = WholeSkyPlotTool
config.atools.psfCModelScatter_psf_cModel_diff_sigmaMad.metric = (
    "psfCModelScatter_{band}_psf_cModel_diff_sigmaMad"
)
config.atools.psfCModelScatter_psf_cModel_diff_mean = WholeSkyPlotTool
config.atools.psfCModelScatter_psf_cModel_diff_mean.metric = (
    "psfCModelScatter_{band}_psf_cModel_diff_mean"
)
config.atools.shapeSizeFractionalDiff_highSNStars_median = WholeSkyPlotTool
config.atools.shapeSizeFractionalDiff_highSNStars_median.metric = (
    "shapeSizeFractionalDiff_{band}_highSNStars_median"
)
config.atools.shapeSizeFractionalDiff_highSNStars_sigmaMad = WholeSkyPlotTool
config.atools.shapeSizeFractionalDiff_highSNStars_sigmaMad.metric = (
    "shapeSizeFractionalDiff_{band}_highSNStars_sigmaMad"
)
config.atools.shapeSizeFractionalDiff_lowSNStars_median = WholeSkyPlotTool
config.atools.shapeSizeFractionalDiff_lowSNStars_median.metric = (
    "shapeSizeFractionalDiff_{band}_lowSNStars_median"
)
config.atools.shapeSizeFractionalDiff_lowSNStars_sigmaMad = WholeSkyPlotTool
config.atools.shapeSizeFractionalDiff_lowSNStars_sigmaMad.metric = (
    "shapeSizeFractionalDiff_{band}_lowSNStars_sigmaMad"
)
config.atools.skyFluxStatisticMetric_medianSky = WholeSkyPlotTool
config.atools.skyFluxStatisticMetric_medianSky.metric = "skyFluxStatisticMetric_{band}_medianSky"
config.atools.skyFluxStatisticMetric_meanSky = WholeSkyPlotTool
config.atools.skyFluxStatisticMetric_meanSky.metric = "skyFluxStatisticMetric_{band}_meanSky"
config.atools.skyFluxStatisticMetric_stdevSky = WholeSkyPlotTool
config.atools.skyFluxStatisticMetric_stdevSky.metric = "skyFluxStatisticMetric_{band}_stdevSky"
config.atools.skyFluxStatisticMetric_sigmaMADSky = WholeSkyPlotTool
config.atools.skyFluxStatisticMetric_sigmaMADSky.metric = "skyFluxStatisticMetric_{band}_sigmaMADSky"

# Keys without `Band``
config.atools.wPerpPSF_wPerp_psfFlux_median = WholeSkyPlotTool
config.atools.wPerpPSF_wPerp_psfFlux_median.metric = "wPerpPSF_wPerp_psfFlux_median"
config.atools.wPerpPSF_wPerp_psfFlux_median.parameterizedBand = False
config.atools.wPerpPSF_wPerp_psfFlux_sigmaMAD = WholeSkyPlotTool
config.atools.wPerpPSF_wPerp_psfFlux_sigmaMAD.metric = "wPerpPSF_wPerp_psfFlux_sigmaMAD"
config.atools.wPerpPSF_wPerp_psfFlux_sigmaMAD.parameterizedBand = False
config.atools.yPerpPSF_yPerp_psfFlux_median = WholeSkyPlotTool
config.atools.yPerpPSF_yPerp_psfFlux_median.metric = "yPerpPSF_yPerp_psfFlux_median"
config.atools.yPerpPSF_yPerp_psfFlux_median.parameterizedBand = False
config.atools.yPerpPSF_yPerp_psfFlux_sigmaMAD = WholeSkyPlotTool
config.atools.yPerpPSF_yPerp_psfFlux_sigmaMAD.metric = "yPerpPSF_yPerp_psfFlux_sigmaMAD"
config.atools.yPerpPSF_yPerp_psfFlux_sigmaMAD.parameterizedBand = False
config.atools.skippedDeblenderMetrics_numSkippedPeaks = WholeSkyPlotTool
config.atools.skippedDeblenderMetrics_numSkippedPeaks.metric = "skippedDeblenderMetrics_numSkippedPeaks"
config.atools.skippedDeblenderMetrics_numSkippedPeaks.parameterizedBand = False
config.atools.skippedDeblenderMetrics_numSkippedBlends = WholeSkyPlotTool
config.atools.skippedDeblenderMetrics_numSkippedBlends.metric = "skippedDeblenderMetrics_numSkippedBlends"
config.atools.skippedDeblenderMetrics_numSkippedBlends.parameterizedBand = False
config.atools.skippedDeblenderMetrics_numBlendParentTooBig = WholeSkyPlotTool
config.atools.skippedDeblenderMetrics_numBlendParentTooBig.metric = (
    "skippedDeblenderMetrics_numBlendParentTooBig"
)
config.atools.skippedDeblenderMetrics_numBlendParentTooBig.parameterizedBand = False
config.atools.skippedDeblenderMetrics_numBlendTooManyPeaks = WholeSkyPlotTool
config.atools.skippedDeblenderMetrics_numBlendTooManyPeaks.metric = (
    "skippedDeblenderMetrics_numBlendTooManyPeaks"
)
config.atools.skippedDeblenderMetrics_numBlendTooManyPeaks.parameterizedBand = False
config.atools.skippedDeblenderMetrics_numBlendTooManyMasked = WholeSkyPlotTool
config.atools.skippedDeblenderMetrics_numBlendTooManyMasked.metric = (
    "skippedDeblenderMetrics_numBlendTooManyMasked"
)
config.atools.skippedDeblenderMetrics_numBlendTooManyMasked.parameterizedBand = False
