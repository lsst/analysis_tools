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

import lsst.utils.tests
import numpy as np
from lsst.analysis.tools.actions.plot.plotUtils import perpDistance, shorten_list, stellarLocusFit


class FitTest(lsst.utils.tests.TestCase):
    """Test to see if the fitting and distance calculations are working."""

    def testFitLine(self):
        """Make a known array of points for x and y and then test that
        the derived fit parameters are as expected."""

        xs = np.arange(1, 10)
        ys = np.arange(1, 10)

        # Define an initial fit box that encompasses the points
        testParams = {"xMin": 0, "xMax": 11, "yMin": 0, "yMax": 11, "mHW": 1, "bHW": 0}
        paramsOut = stellarLocusFit(xs, ys, testParams)

        # stellarLocusFit performs two iterations of fitting and also
        # calculates the perpendicular gradient to the fit line and
        # the points of intersection between the box and the fit
        # line. Test that these are returning what is expected.
        self.assertFloatsAlmostEqual(paramsOut["mODR"], 1.0)
        self.assertFloatsAlmostEqual(paramsOut["mODR2"], 1.0)
        self.assertFloatsAlmostEqual(paramsOut["bODR"], 0.0)
        self.assertFloatsAlmostEqual(paramsOut["bODR2"], 0.0)
        self.assertFloatsAlmostEqual(paramsOut["bPerpMin"], 0.0)
        self.assertFloatsAlmostEqual(paramsOut["bPerpMax"], 22.0)

    def testPerpDistance(self):
        """Test the calculation of the perpendicular distance."""

        p1 = np.array([1, 1])
        p2 = np.array([2, 2])
        testPoints = np.array([[1, 2], [1.5, 1.5], [2, 1]])
        # perpDistance uses two points, p1 and p2, to define a line
        # then calculates the perpendicular distance of the testPoints
        # to this line
        dists = perpDistance(p1, p2, testPoints)
        self.assertFloatsAlmostEqual(np.array([1.0 / np.sqrt(2), 0.0, -1.0 / np.sqrt(2)]), np.array(dists))


class ShortListTestCase(lsst.utils.tests.TestCase):
    """Test to see if the shorten_list function works as it should."""

    def test_shorten_list(self):
        self.assertEqual(shorten_list([]), "")  # empty container
        self.assertEqual(shorten_list([5]), "5")  # single element
        self.assertEqual(shorten_list([5, 6]), "5-6")  # 2 contigous elements
        self.assertEqual(shorten_list([5, 7]), "5,7")  # 2 non-contiguous elements
        self.assertEqual(shorten_list([-7, -6, -5]), "-7--5")  # 3 contiguous elements
        self.assertEqual(shorten_list([5, 7, 9]), "5,7,9")  # 3 non-contiguous elements
        self.assertEqual(shorten_list([5, 6, 8]), "5-6,8")  # 3 mixed elements

        # Test for different data types
        self.assertEqual(shorten_list({1, 2, 3, 5, 7, 8, 9, 10, 13}), "1-3,5,7-10,13")  # test for a set
        self.assertEqual(
            shorten_list((1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18)), "1-3,5-13,15-18"
        )  # test for a tuple
        self.assertEqual(shorten_list(range(100)), "0-99")  # test for an itertor
        self.assertEqual(shorten_list(np.arange(100, dtype=int)), "0-99")  # test for a numpy array

        # Test the RC2/DC2 tract list explicitly.
        self.assertEqual(shorten_list([9615, 9813, 9697]), "9615,9697,9813")
        self.assertEqual(shorten_list([3828, 3829]), "3828-3829")

        # Test the keyword-only arguments.
        self.assertEqual(shorten_list([1, 2, 3, 4, 7, 8, 9], range_indicator=".."), "1..4,7..9")
        self.assertEqual(shorten_list([1, 2, 3, 4, 7, 8, 9], range_separator="^"), "1-4^7-9")
        self.assertEqual(
            shorten_list([1, 2, 3, 4, 7, 8, 9], range_indicator=":", range_separator=";"), "1:4;7:9"
        )

        # Test that those are keyword only.
        with self.assertRaises(TypeError):
            shorten_list([1, 2, 3, 4, 7, 8, 9], "..", "^")
        with self.assertRaises(TypeError):
            shorten_list([1, 2, 3, 4, 7, 8, 9], "..", range_separator="^")
        with self.assertRaises(TypeError):
            shorten_list([1, 2, 3, 4, 7, 8, 9], "!", range_indicator="..")

        # Test that it sorts the container.
        self.assertEqual(shorten_list([6, 3, 1, 2]), "1-3,6")

        # Repeated numbers are not expected. Test that the function can
        # nevertheless handle it.
        self.assertEqual(shorten_list([1, 2, 2, 3, 3, 3, 5, 7, 8, 7]), "1-3,5,7-8")


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
