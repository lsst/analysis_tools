# Configuration for instances of lsst.analysis.tools.tasks.WholeSkyAnalysisTask.
# Generated whole sky plots showsing per-tract aggregations of coadd mask
# fractions.

from lsst.analysis.tools.atools import WholeSkyPlotTool

maskPlanes = {"NO_DATA", "CLIPPED", "NOT_DEBLENDED", "DETECTED"}

if hasattr(parameters, "exclude_aggregate_mask_planes"):
    maskPlanes -= set(parameters.exclude_aggregate_mask_planes)

for maskPlane in maskPlanes:
    atoolName = f"{maskPlane}_fraction_mean"
    setattr(config.atools, atoolName, WholeSkyPlotTool)
    atool = getattr(config.atools, atoolName)

    if maskPlane == "NO_DATA":
        setattr(atool, "metric", f"full_tract_mean_{maskPlane}_fraction")
    else:
        setattr(atool, "metric", f"full_tract_mean_{maskPlane}_valid_data_fraction")

    setattr(atool, "parameterizedBand", False)
    setattr(atool.produce.plot, "showOutliers", False)
    setattr(atool.produce.plot, "colorMapType", "sequential")
