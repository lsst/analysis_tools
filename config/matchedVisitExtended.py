# Configuration for extended instances of
# lsst.analysis.tools.tasks.AssociatedSourcesTractAnalysisTask
#
# This does not include the Core configuration, as it may be used by instances
# that are configured to optionally run alongside instances with that
# configurations.

from lsst.analysis.tools.atools import *
from lsst.analysis.tools.interfaces import *

config.atools.modelPhotRepStarSn5to10 = StellarPhotometricRepeatability
config.atools.modelPhotRepStarSn5to10.fluxType = "gaussianFlux"
config.atools.modelPhotRepStarSn5to10.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 5
config.atools.modelPhotRepStarSn5to10.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 10
config.atools.modelPhotRepStarSn5to10.produce.plot = NoPlot

config.atools.modelPhotRepStarSn10to20 = StellarPhotometricRepeatability
config.atools.modelPhotRepStarSn10to20.fluxType = "gaussianFlux"
config.atools.modelPhotRepStarSn10to20.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 10
config.atools.modelPhotRepStarSn10to20.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 20
config.atools.modelPhotRepStarSn10to20.produce.plot = NoPlot

config.atools.modelPhotRepStarSn20to40 = StellarPhotometricRepeatability
config.atools.modelPhotRepStarSn20to40.fluxType = "gaussianFlux"
config.atools.modelPhotRepStarSn20to40.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 20
config.atools.modelPhotRepStarSn20to40.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 40
config.atools.modelPhotRepStarSn20to40.produce.plot = NoPlot

config.atools.modelPhotRepStarSn40to80 = StellarPhotometricRepeatability
config.atools.modelPhotRepStarSn40to80.fluxType = "gaussianFlux"
config.atools.modelPhotRepStarSn40to80.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 40
config.atools.modelPhotRepStarSn40to80.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 80
config.atools.modelPhotRepStarSn40to80.produce.plot = NoPlot

config.atools.psfPhotRepStarSn5to10 = StellarPhotometricRepeatability
config.atools.psfPhotRepStarSn5to10.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 5
config.atools.psfPhotRepStarSn5to10.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 10
config.atools.psfPhotRepStarSn5to10.produce.plot = NoPlot

config.atools.psfPhotRepStarSn10to20 = StellarPhotometricRepeatability
config.atools.psfPhotRepStarSn10to20.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 10
config.atools.psfPhotRepStarSn10to20.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 20
config.atools.psfPhotRepStarSn10to20.produce.plot = NoPlot

config.atools.psfPhotRepStarSn20to40 = StellarPhotometricRepeatability
config.atools.psfPhotRepStarSn20to40.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 20
config.atools.psfPhotRepStarSn20to40.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 40
config.atools.psfPhotRepStarSn20to40.produce.plot = NoPlot

config.atools.psfPhotRepStarSn40to80 = StellarPhotometricRepeatability
config.atools.psfPhotRepStarSn40to80.process.filterActions.perGroupStdevFiltered.selectors.sn.minimum = 40
config.atools.psfPhotRepStarSn40to80.process.filterActions.perGroupStdevFiltered.selectors.sn.maximum = 80
config.atools.psfPhotRepStarSn40to80.produce.plot = NoPlot
