# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
from __future__ import annotations

__all__ = (
    "CoaddInputAnalysisConfig",
    "CoaddInputAnalysisTask",
)

from collections import defaultdict

import astropy.table
import numpy as np

from lsst.geom import SpherePoint, degrees
from lsst.pipe.base import (
    InputQuantizedConnection,
    OutputQuantizedConnection,
    QuantumContext,
)
from lsst.pipe.base import connectionTypes as cT
from lsst.skymap import BaseSkyMap
from lsst.sphgeom import ConvexPolygon

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CoaddInputAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("tract", "patch", "band"),
    defaultTemplates={"outputName": "coaddInputAnalysis"},
):
    raw = cT.Input(
        doc="Raw exposure data.",
        name="raw",
        storageClass="Exposure",
        dimensions=(
            "instrument",
            "exposure",
            "detector",
        ),
        deferLoad=True,
        multiple=True,
    )
    coaddInputs = cT.Input(
        doc="List of dataIds that went into a coadd.",
        name="deep_coadd_input_summary_tract",
        storageClass="ArrowAstropy",
        dimensions=(
            "skymap",
            "tract",
            "band",
        ),
    )
    visitTable = cT.Input(
        doc="""Table summarising metadata at detector-level for all visist,
            from which the per-detector image corner coordinates are used
            to determine which calibrated images spatially overlap a given
            patch.""",
        name="ccdVisitTable",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
        multiple=False,
    )
    visitSummary = cT.Input(
        doc="""Per-visit summary of calibrated images. Not directly used in
            processing, but including it as a per-visit input substantially
            reduces the time to build the quantum graph.""",
        name="finalVisitSummary",
        storageClass="ExposureCatalog",
        dimensions=(
            "instrument",
            "visit",
        ),
        deferLoad=True,
        multiple=True,
    )
    skyMap = cT.PrerequisiteInput(
        doc="Sky map defining the tracts and patches.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )


class CoaddInputAnalysisConfig(AnalysisBaseConfig, pipelineConnections=CoaddInputAnalysisConnections):
    pass


class CoaddInputAnalysisTask(AnalysisPipelineTask):
    """
    Construct a table containing the visit and detector IDs
    of all PVIs that the butler believes could potentially
    have gone into a coadd. Indicate those whose post-calibration
    astrometry overlaps with the patch. Also Indicate those which
    finally made it into the coadd.
    """

    ConfigClass = CoaddInputAnalysisConfig
    _DefaultName = "coaddInputAnalysis"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:

        # The PVIs that actually made it into the coadd:
        coaddInputTable = butlerQC.get(inputRefs.coaddInputs)
        patch = butlerQC.quantum.dataId["patch"]
        patchInputs = coaddInputTable[coaddInputTable["patch"] == patch]
        inCoadd = set(zip(patchInputs["visit"], patchInputs["detector"]))

        # On-sky polygon of patch to determine which PVIs cover it.
        dataId = butlerQC.quantum.dataId
        skyMap = butlerQC.get(inputRefs.skyMap)
        tractInfo = skyMap[dataId["tract"]]
        patchInfo = tractInfo.getPatchInfo(dataId["patch"])
        patchPoly = patchInfo.getOuterSkyPolygon()

        visitTableHandle = butlerQC.get(inputRefs.visitTable)
        columns = [
            "visitId",
            "detector",
            "llcra",
            "llcdec",
            "ulcra",
            "ulcdec",
            "urcra",
            "urcdec",
            "lrcra",
            "lrcdec",
        ]
        imageCorners = self.loadData(visitTableHandle, columns)

        # Group raws by visit so the imageCorners table can be reduced by visit
        rawsByVisit = defaultdict(list)
        for rawRef in inputRefs.raw:
            rawsByVisit[rawRef.dataId["exposure"]].append(rawRef)

        data = self.makeData(inCoadd, patchPoly, imageCorners, rawsByVisit)

        outputs = self.run(data=data)
        butlerQC.put(outputs, outputRefs)

    def makeData(self, inCoadd, patchPoly, imageCorners, rawsByVisit):
        """Build the per-PVI data table for analysis.

        Parameters
        ----------
        inCoadd : `set` [`tuple` [`int`, `int`]]
            Set of (visit, detector) pairs that made it into the coadd.
        patchPoly : `lsst.sphgeom.ConvexPolygon`
            On-sky polygon of the patch.
        imageCorners : `astropy.table.Table`
            Table of per-detector image corner coordinates.
        rawsByVisit : `dict` [`int`, `list`]
            Raw data references grouped by visit ID.

        Returns
        -------
        data : `astropy.table.Table`
            Table with columns ``visit``, ``detector``,
            ``visitSummaryRecord``, ``patchOverlap``, and ``inCoadd``.
        """
        n = sum(len(rawRefs) for rawRefs in rawsByVisit.values())
        visits = [None] * n
        detectors = [None] * n
        inCoadd_col = [False] * n
        visitSummaryRecord = [False] * n
        overlapsWithPatch = [False] * n

        i = 0
        for visit, rawRefs in rawsByVisit.items():
            visitSummary = imageCorners[imageCorners["visitId"] == visit]

            if len(visitSummary) == 0:
                for rawRef in rawRefs:
                    visits[i] = visit
                    detectors[i] = rawRef.dataId["detector"]
                    i += 1
                continue

            for rawRef in rawRefs:
                detector = rawRef.dataId["detector"]
                visits[i] = visit
                detectors[i] = detector
                inCoadd_col[i] = (visit, detector) in inCoadd

                if not np.any(rowMask := visitSummary["detector"] == detector):
                    # If there's no visit summary, move on.
                    # visitSummaryRecord and overlapsWithPatch already contain
                    # False as default:
                    i += 1
                    continue

                # If we have got this far, then there is a row for this image.
                # Record whether the calibrated image covers the patch:
                visitSummaryRecord[i] = True

                # It is a Bad Thing if there is more than one row per
                # (visit,detector):
                if (nMatches := rowMask.sum()) > 1:
                    raise RuntimeError(
                        f"Expected exactly one row for visit={visit}, "
                        f"detector={detector} in imageCorners, found {nMatches}."
                    )
                raCorners = np.array(
                    [
                        visitSummary[rowMask][corner].value[0]
                        for corner in ["llcra", "ulcra", "urcra", "lrcra"]
                    ]
                )
                decCorners = np.array(
                    [
                        visitSummary[rowMask][corner].value[0]
                        for corner in ["llcdec", "ulcdec", "urcdec", "lrcdec"]
                    ]
                )

                # Defensively catch non-finite corner coords. Treat such
                # detectors conservatively as not overlapping the patch.
                if np.all(np.isfinite(raCorners)) and np.all(np.isfinite(decCorners)):
                    detCorners = [
                        SpherePoint(ra, dec, units=degrees).getVector()
                        for ra, dec in zip(raCorners, decCorners)
                    ]
                    detPoly = ConvexPolygon.convexHull(detCorners)
                    overlapsWithPatch[i] = patchPoly.intersects(detPoly)

                i += 1

        return astropy.table.Table(
            {
                "visit": visits,
                "detector": detectors,
                "visitSummaryRecord": visitSummaryRecord,
                "patchOverlap": overlapsWithPatch,
                "inCoadd": inCoadd_col,
            }
        )
