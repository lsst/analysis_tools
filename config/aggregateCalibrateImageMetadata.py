# Configuration for lsst.analysis.tools.tasks.AggregatedTaskMetadataAnalysisTask
# when aggregating calibrateImage task metadata across detectors to produce
# per-visit aggregate statistics.

from lsst.analysis.tools.actions.scalar import CountAction, MaxAction, MedianAction, MinAction, SigmaMadAction
from lsst.analysis.tools.atools import AggregatedTaskMetadataMetricTool

config.atools.calibrateImageMetadataMetrics = AggregatedTaskMetadataMetricTool
config.atools.calibrateImageMetadataMetrics.taskName = "calibrateImage"
config.atools.calibrateImageMetadataMetrics.metrics = {
    "initial_psf_positive_footprint_count": "ct",
    "initial_psf_negative_footprint_count": "ct",
    "initial_psf_positive_peak_count": "ct",
    "initial_psf_negative_peak_count": "ct",
    "simple_psf_positive_footprint_count": "ct",
    "simple_psf_negative_footprint_count": "ct",
    "simple_psf_positive_peak_count": "ct",
    "simple_psf_negative_peak_count": "ct",
    "bad_mask_fraction": "",
    "cr_mask_fraction": "",
    "crosstalk_mask_fraction": "",
    "detected_mask_fraction": "",
    "detected_negative_mask_fraction": "",
    "edge_mask_fraction": "",
    "intrp_mask_fraction": "",
    "no_data_mask_fraction": "",
    "sat_mask_fraction": "",
    "suspect_mask_fraction": "",
    "unmaskednan_mask_fraction": "",
    "numAvailStars": "ct",
    "numGoodStars": "ct",
    "sky_footprint_count": "ct",
    "post_deblend_source_count": "ct",
    "star_count": "ct",
    "saturated_source_count": "ct",
    "bad_source_count": "ct",
    "cosmic_ray_count": "ct",
    "matched_psf_star_count": "ct",
    "final_psf_sigma": "pixel",
    "astrometry_matches_count": "ct",
    "photometry_matches_count": "ct",
    "bg_subtracted_skyPixel_instFlux_median": "adu",
    "bg_subtracted_skyPixel_instFlux_stdev": "adu",
    "bg_subtracted_skySource_flux_median": "nJy",
    "bg_subtracted_skySource_flux_stdev": "nJy",
    "adaptive_threshold_value": "",
    "initial_to_final_wcs": "arcsec",
    "astrom_offset_mean": "arcsec",
    "astrom_offset_std": "arcsec",
    "astrom_offset_median": "arcsec",
    "failed_deblend_source_count": "ct",
}
config.atools.calibrateImageMetadataMetrics.subTaskNames = {
    "numAvailStars": "psf_measure_psf",
    "numGoodStars": "psf_measure_psf",
    "sky_footprint_count": "star_sky_sources",
    "cosmic_ray_count": "psf_repair",
}
config.atools.calibrateImageMetadataMetrics.newNames = {
    "numAvailStars": "psf_available_star_count",
    "numGoodStars": "psf_good_star_count",
    "adaptive_threshold_value": "final_adaptive_threshold_value",
}

aggregations = {"min": MinAction, "max": MaxAction, "median": MedianAction, "mad": SigmaMadAction, "ct": CountAction}
config.atools.calibrateImageMetadataMetrics.aggregationUnits = {"ct": "ct"}

for metric in config.atools.calibrateImageMetadataMetrics.metrics:
    for agg_name, agg_cls in aggregations.items():
        action = agg_cls()
        action.vectorKey = metric
        setattr(
            config.atools.calibrateImageMetadataMetrics.process.calculateActions,
            f"{metric}_{agg_name}",
            action,
        )
