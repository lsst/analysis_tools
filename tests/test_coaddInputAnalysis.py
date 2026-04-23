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

import unittest

import astropy.table
import numpy as np

import lsst.utils.tests
from lsst.analysis.tools.tasks.coaddInputAnalysis import (
    CoaddInputAnalysisConfig,
    CoaddInputAnalysisTask,
)
from lsst.geom import SpherePoint, degrees
from lsst.sphgeom import ConvexPolygon


class MockRawRef:
    """Minimal stand-in for a butler raw data reference."""

    def __init__(self, detector):
        self.dataId = {"detector": detector}


def makePatchPoly(ra_center, dec_center, half_size=0.5):
    """Return a square ConvexPolygon centred at (ra_center, dec_center)."""
    corners = [
        SpherePoint(ra_center - half_size, dec_center - half_size, degrees).getVector(),
        SpherePoint(ra_center + half_size, dec_center - half_size, degrees).getVector(),
        SpherePoint(ra_center + half_size, dec_center + half_size, degrees).getVector(),
        SpherePoint(ra_center - half_size, dec_center + half_size, degrees).getVector(),
    ]
    return ConvexPolygon.convexHull(corners)


def makeImageCornersRow(visit, detector, ra_center, dec_center, half_size=0.1):
    """Return a single-row imageCorners Table for one detector."""
    return astropy.table.Table(
        {
            "visitId": [visit],
            "detector": [detector],
            "llcra": [ra_center - half_size],
            "llcdec": [dec_center - half_size],
            "ulcra": [ra_center - half_size],
            "ulcdec": [dec_center + half_size],
            "urcra": [ra_center + half_size],
            "urcdec": [dec_center + half_size],
            "lrcra": [ra_center + half_size],
            "lrcdec": [dec_center - half_size],
        }
    )


class TestMakeData(lsst.utils.tests.TestCase):
    """Tests for CoaddInputAnalysisTask.makeData."""

    def setUp(self):
        self.task = CoaddInputAnalysisTask(config=CoaddInputAnalysisConfig())
        # 1-degree square patch centred at (180, 0).
        self.patchPoly = makePatchPoly(180.0, 0.0, half_size=0.5)

    def test_overlap_and_in_coadd(self):
        """PVI inside patch and in the coadd is correctly flagged."""
        visit, detector = 100, 5
        imageCorners = makeImageCornersRow(visit, detector, 180.0, 0.0)
        rawsByVisit = {visit: [MockRawRef(detector)]}
        inCoadd = {(visit, detector)}

        data = self.task.makeData(inCoadd, self.patchPoly, imageCorners, rawsByVisit)

        self.assertEqual(len(data), 1)
        self.assertTrue(data["patchOverlap"][0])
        self.assertTrue(data["visitSummaryRecord"][0])
        self.assertTrue(data["inCoadd"][0])

    def test_no_overlap_not_in_coadd(self):
        """PVI outside patch and not in the coadd is correctly flagged."""
        visit, detector = 100, 5
        imageCorners = makeImageCornersRow(visit, detector, 190.0, 0.0)
        rawsByVisit = {visit: [MockRawRef(detector)]}

        data = self.task.makeData(set(), self.patchPoly, imageCorners, rawsByVisit)

        self.assertEqual(len(data), 1)
        self.assertFalse(data["patchOverlap"][0])
        self.assertTrue(data["visitSummaryRecord"][0])
        self.assertFalse(data["inCoadd"][0])

    def test_non_finite_corners_treated_as_no_overlap(self):
        """
        A PVI with a non-finite corner coordinate is treated as
        not overlapping.
        """
        visit, detector = 100, 5
        imageCorners = makeImageCornersRow(visit, detector, 180.0, 0.0)
        imageCorners["llcra"][0] = np.nan
        rawsByVisit = {visit: [MockRawRef(detector)]}

        data = self.task.makeData(set(), self.patchPoly, imageCorners, rawsByVisit)

        self.assertFalse(data["patchOverlap"][0])
        self.assertTrue(data["visitSummaryRecord"][0])

    def test_visit_absent_from_image_corners(self):
        """A visit with no imageCorners entry gets all-False flags."""
        visit, detector = 100, 5
        imageCorners = astropy.table.Table(
            {
                col: []
                for col in [
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
            }
        )
        rawsByVisit = {visit: [MockRawRef(detector)]}

        data = self.task.makeData(set(), self.patchPoly, imageCorners, rawsByVisit)

        self.assertEqual(len(data), 1)
        self.assertFalse(data["visitSummaryRecord"][0])
        self.assertFalse(data["patchOverlap"][0])
        self.assertFalse(data["inCoadd"][0])

    def test_detector_absent_from_visit_corners(self):
        """
        A detector with no row in imageCorners gets False for corner-derived
        flags.
        """
        visit, detector = 100, 5
        other_detector = 6
        imageCorners = makeImageCornersRow(visit, other_detector, 180.0, 0.0)
        rawsByVisit = {visit: [MockRawRef(detector)]}

        data = self.task.makeData(set(), self.patchPoly, imageCorners, rawsByVisit)

        self.assertFalse(data["visitSummaryRecord"][0])
        self.assertFalse(data["patchOverlap"][0])

    def test_duplicate_detector_row_raises(self):
        """
        Duplicate (visit, detector) rows in imageCorners raise RuntimeError.
        """
        visit, detector = 100, 5
        row = makeImageCornersRow(visit, detector, 180.0, 0.0)
        imageCorners = astropy.table.vstack([row, row])
        rawsByVisit = {visit: [MockRawRef(detector)]}

        with self.assertRaises(RuntimeError):
            self.task.makeData(set(), self.patchPoly, imageCorners, rawsByVisit)

    def test_multiple_visits_and_detectors(self):
        """
        Multiple visits and detectors produce the correct number of rows with
        correct flags.
        """
        v1, v2 = 100, 200
        d1, d2, d3 = 5, 6, 7
        imageCorners = astropy.table.vstack(
            [
                makeImageCornersRow(v1, d1, 180.0, 0.0),  # overlaps
                makeImageCornersRow(v1, d2, 190.0, 0.0),  # no overlap
                makeImageCornersRow(v2, d3, 180.0, 0.0),  # overlaps
            ]
        )
        rawsByVisit = {
            v1: [MockRawRef(d1), MockRawRef(d2)],
            v2: [MockRawRef(d3)],
        }
        inCoadd = {(v1, d1), (v2, d3)}

        data = self.task.makeData(inCoadd, self.patchPoly, imageCorners, rawsByVisit)
        data.sort(["visit", "detector"])

        self.assertEqual(len(data), 3)

        row_v1d1 = data[(data["visit"] == v1) & (data["detector"] == d1)][0]
        self.assertTrue(row_v1d1["patchOverlap"])
        self.assertTrue(row_v1d1["inCoadd"])

        row_v1d2 = data[(data["visit"] == v1) & (data["detector"] == d2)][0]
        self.assertFalse(row_v1d2["patchOverlap"])
        self.assertFalse(row_v1d2["inCoadd"])

        row_v2d3 = data[(data["visit"] == v2) & (data["detector"] == d3)][0]
        self.assertTrue(row_v2d3["patchOverlap"])
        self.assertTrue(row_v2d3["inCoadd"])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
