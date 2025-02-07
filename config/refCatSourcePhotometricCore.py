# Configuration for standard instances of
# lsst.analysis.tools.tasks.refCatSourceAnalysis.RefCatSourceAnalysisTask
#
# These configurations are not in the class default configs because there may
# be instances of the class that are intended to optionally run alongside one
# that includes this standard set.

from lsst.analysis.tools.atools import *

config.atools.photomDiffPsfSkyVisitPlot = TargetRefCatDeltaPsfSkyVisitPlot
config.atools.photomDiffAp09SkyVisitPlot = TargetRefCatDeltaAp09SkyVisitPlot
config.atools.photoDiffPsfScatterVisitPlot = TargetRefCatDeltaPsfScatterVisitPlot
config.atools.photoDiffCModelScatterVisitPlot = TargetRefCatDeltaAp09ScatterVisitPlot
