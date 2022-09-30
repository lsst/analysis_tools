from lsst.analysis.tools.actions.plot.histPlot import HistPlot, HistPanel
from lsst.analysis.tools.interfaces import AnalysisPlot
from ..analysisParts.demoPsfApRatio import DemoPsfApRatioBaseClass


class DemoPlot(AnalysisPlot, DemoPsfApRatioBaseClass):
    def setDefaults(self):
        super().setDefaults()

        self.produce = HistPlot()

        self.produce.panels["panel_flux"] = HistPanel()
        self.produce.panels["panel_flux"].label = "Psf/Ap Ratio"
        self.produce.panels["panel_flux"].hists = dict(fluxRatioMetric="fluxRatioMetric")
