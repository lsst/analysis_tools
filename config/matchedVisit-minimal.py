# Configuration for minimal instances of
# lsst.analysis.tools.tasks.AssociatedSourcesTractAnalysisTask
#
# These configurations are not in the class default configs because there may
# be instances of the class that are intended to optionally run alongside one
# that includes this minimal set.

from lsst.analysis.tools.atools import *

config.atools.stellarAstrometricSelfRepeatabilityRA = AstrometricRepeatability
config.atools.stellarAstrometricSelfRepeatabilityRA.coordinate = "RA"
config.atools.stellarAstrometricSelfRepeatabilityDec = AstrometricRepeatability
config.atools.stellarAstrometricSelfRepeatabilityDec.coordinate = "Dec"
