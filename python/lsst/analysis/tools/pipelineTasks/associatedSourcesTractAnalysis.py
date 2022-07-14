from __future__ import annotations


from lsst.pipe.base import connectionTypes as ct

from .base import AnalysisBaseConnections, AnalysisBaseConfig, AnalysisPipelineTask

# These need to be updated for this analysis context
#from ..analysisPlots.analysisPlots import ShapeSizeFractionalDiffScatter, Ap12_PSF_skyPlot
#from ..analysisMetrics.analysisMetrics import ShapeSizeFractionalMetric


class AssociatedSourcesTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract"),
    defaultTemplates={"associatedSourcesInputName": "isolated_star_sources"},
):
    data = ct.Input(
        doc="Visit based source table to load from the butler",
        name="sourceTable_visit",
        storageClass="DataFrame",
        #deferLoad=True,
        dimensions=("visit", "band"),
        mutliple=True,
    )

    associated_sources = ct.Input(
        doc="Table of associated sources",
        name="{associatedSourcesInputName}",
        storageClass="DataFrame",
        #deferLoad=True,
        dimentions=("skymap", "tract"),
    )

class AssociatedSourcesTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=AssociatedSourcesTractAnalysisConnections
):
    def setDefaults(self):
        super().setDefaults()
        # set plots to run
        # update for this analysis context
        #self.plots.shapeSizeFractionalDiffScatter = ShapeSizeFractionalDiffScatter()
        #self.plots.Ap12_PSF_skyPlot = Ap12_PSF_skyPlot()

        # set metrics to run
        # update for this analysis context
        #self.metrics.shapeSizeFractionalMetric = ShapeSizeFractionalMetric()


class AssociatedSourcesTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = AssociatedSourcesTractAnalysisConfig
    _DefaultName = "associatedSourcesTractAnalysisTask"