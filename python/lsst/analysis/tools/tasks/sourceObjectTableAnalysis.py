from __future__ import annotations

__all__ = ("SourceObjectTableAnalysisConfig", "SourceObjectTableAnalysisTask")

import lsst.pex.config as pexConfig
import numpy as np
import pandas as pd
from astropy.table import vstack
from lsst.pipe.base import connectionTypes as ct
from smatch import Matcher

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class SourceObjectTableAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit",),
    defaultTemplates={
        "inputName": "sourceTable_visit",
        "inputCoaddName": "deep",
        "associatedSourcesInputName": "isolated_star_sources",
        "outputName": "sourceObjectTable",
    },
):
    data = ct.Input(
        doc="Visit based source table to load from the butler",
        name="sourceTable_visit",
        storageClass="ArrowAstropy",
        dimensions=("visit", "band"),
        deferLoad=True,
    )

    associatedSources = ct.Input(
        doc="Table of associated sources",
        name="{associatedSourcesInputName}",
        storageClass="ArrowAstropy",
        multiple=True,
        deferLoad=True,
        dimensions=("instrument", "skymap", "tract"),
    )

    refCat = ct.Input(
        doc="Catalog of positions to use as reference.",
        name="objectTable",
        storageClass="DataFrame",
        dimensions=["skymap", "tract", "patch"],
        multiple=True,
        deferLoad=True,
    )


class SourceObjectTableAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=SourceObjectTableAnalysisConnections
):
    ra_column = pexConfig.Field(
        doc="Name of column in refCat to use for right ascension.",
        dtype=str,
        default="r_ra",
    )
    dec_column = pexConfig.Field(
        doc="Name of column in refCat to use for declination.",
        dtype=str,
        default="r_dec",
    )


class SourceObjectTableAnalysisTask(AnalysisPipelineTask):
    ConfigClass = SourceObjectTableAnalysisConfig
    _DefaultName = "sourceTableVisitAnalysis"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        # Get isolated sources:
        visit = inputs["data"].dataId["visit"]
        band = inputs["data"].dataId["band"]
        names = self.collectInputNames()
        names -= {self.config.ra_column, self.config.dec_column}
        data = inputs["data"].get(parameters={"columns": names})

        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)

        isolatedSources = []
        for associatedSourcesRef in inputs["associatedSources"]:
            associatedSources = associatedSourcesRef.get(parameters={"columns": ["visit", "source_row"]})
            visit_sources = associatedSources[associatedSources["visit"] == visit]
            isolatedSources.append(data[visit_sources["source_row"]])
        isolatedSources = vstack(isolatedSources)

        # Get objects:
        allRefCats = []
        for refCatRef in inputs["refCat"]:
            refCat = refCatRef.get(
                parameters={"columns": ["detect_isPrimary", self.config.ra_column, self.config.dec_column]}
            )
            goodInds = (
                refCat["detect_isPrimary"]
                & np.isfinite(refCat[self.config.ra_column])
                & np.isfinite(refCat[self.config.dec_column])
            )
            allRefCats.append(refCat[goodInds])

        refCat = pd.concat(allRefCats)

        with Matcher(isolatedSources["coord_ra"], isolatedSources["coord_dec"]) as m:
            idx, isolatedMatchIndices, refMatchIndices, dists = m.query_radius(
                refCat[self.config.ra_column].values,
                refCat[self.config.dec_column].values,
                1 / 3600.0,
                return_indices=True,
            )

        matchRef = refCat.iloc[refMatchIndices]
        matchIS = isolatedSources[isolatedMatchIndices].to_pandas()

        allCat = pd.concat([matchRef.reset_index(), matchIS.reset_index()], axis=1)
        outputs = self.run(data=allCat, bands=band, plotInfo=plotInfo)
        butlerQC.put(outputs, outputRefs)
