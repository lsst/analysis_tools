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

__all__ = ("VisitDetectorTractAnalysisConfig", "VisitDetectorTractAnalysisTask")

import numpy as np
from lsst.geom import Box2D
from lsst.pipe.base import connectionTypes as ct
from lsst.pipe.base import NoWorkFound
from lsst.skymap import BaseSkyMap
from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class VisitDetectorTractAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap", "tract", "instrument"),
    defaultTemplates={
        "outputName": "visit_detector_summary",
    },
):
    visitDetectorTable = ct.Input(
        doc="Table containing detector-level summary statistics for processed visit images.",
        name="preliminary_visit_detector_table",
        storageClass="ArrowAstropy",
        dimensions=("instrument",),
        deferLoad=True,
    )

    skyMap = ct.Input(
        doc="Input definition of geometry/bbox and projection/wcs for warped exposures",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )


class VisitDetectorTractAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=VisitDetectorTractAnalysisConnections
):
    pass


class VisitDetectorTractAnalysisTask(AnalysisPipelineTask):
    ConfigClass = VisitDetectorTractAnalysisConfig
    _DefaultName = "visitDetectorTractAnalysis"

    @staticmethod
    def getBoxWcs(skymap, tract):
        """Get box that defines tract boundaries."""
        tractInfo = skymap.generateTract(tract)
        wcs = tractInfo.getWcs()
        tractBox = tractInfo.getBBox()
        return tractBox, wcs

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)

        dataId = butlerQC.quantum.dataId

        names = self.collectInputNames()
        names |= {'detectorId', 'visitId'}
        cornerPairs = [('llcra', 'llcdec'),
                       ('ulcra', 'ulcdec'),
                       ('urcra', 'urcdec'),
                       ('lrcra', 'lrcdec'),]
        names.update(cornerKey for cornerPair in cornerPairs for cornerKey in cornerPair)

        curatedTable = self.loadData(inputs['visitDetectorTable'], names)

        box, wcs = self.getBoxWcs(inputs['skyMap'], dataId['tract'])
        box = Box2D(box)
        inTract = np.full(len(curatedTable), False)
        for cornerPair in cornerPairs:
            ra = np.radians(curatedTable[cornerPair[0]])
            dec = np.radians(curatedTable[cornerPair[1]])
            x, y = wcs.skyToPixelArray(ra, dec)
            inTract |= box.contains(x, y)

        if np.sum(inTract) == 0:
            raise NoWorkFound(f"No visit detector images overlap tract {dataId.tract.id}")

        data = curatedTable[inTract]

        outputs = self.run(data=data)

        self.putByBand(butlerQC, outputs, outputRefs)

