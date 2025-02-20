# Configuration for standard instances of
# lsst.analysis.tools.tasks.SourceTableVisitAnalysisTask
#
# These configurations are not in the class default configs because there may
# be instances of the class that are intended to optionally run alongside one
# that includes this standard set.

from lsst.analysis.tools.atools import *
from lsst.analysis.tools.contexts import VisitContext

config.atools.skyFluxVisitStatisticMetric = SkyFluxStatisticMetric
config.atools.skyFluxVisitStatisticMetric.applyContext = VisitContext
config.atools.skySourceSky = SkySourceSkyPlot
config.atools.skySourceFlux = SkySourceHistPlot
config.atools.relativeSizeResidualPlot = RelativeSizeResidualPlot
config.addOutputNamePrefix = True
