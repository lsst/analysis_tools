# Configuration for standard instances of
# lsst.analysis.tools.tasks.refCatSourceAnalysis.RefCatSourceAnalysisTask
#
# These configurations are not in the class default configs because there may
# be instances of the class that are intended to optionally run alongside one
# that includes this standard set.

from lsst.analysis.tools.atools import *
from lsst.analysis.tools.contexts import VisitContext

config.atools.astromDiffRASkyVisitPlot = TargetRefCatDeltaRASkyVisitPlot
config.atools.astromDiffDecSkyVisitPlot = TargetRefCatDeltaDecSkyVisitPlot
config.atools.astromDiffRAScatterVisitPlot = TargetRefCatDeltaRAScatterVisitPlot
config.atools.astromDiffDecScatterVisitPlot = TargetRefCatDeltaDecScatterVisitPlot
config.atools.astromDiffMetrics = TargetRefCatDeltaMetrics
config.atools.astromDiffMetrics.applyContext = VisitContext
