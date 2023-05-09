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

import galsim
import lsst.utils.tests
import numpy as np
from lsst.afw.geom import Quadrupole
from lsst.analysis.tools.actions.vector import CalcE, CalcE1, CalcE2, CalcEDiff, CalcMomentSize
from lsst.pex.config import FieldValidationError


class ShapeSizeTestCase(lsst.utils.tests.TestCase):
    """Test ellipiticity and size calculations."""

    @classmethod
    def setUpClass(cls):
        cls.data = np.array(
            [
                (1.3, 1.3, 0.0),  # e1 = e2 = 0
                (2.4, 1.2, 0.6),  # e1 = e2 != 0
                (1.0, 2.0, 0.0),  # e1 < 0; e2 = 0
                (3.5, 3.5, 0.5),  # e1 = 0; e2 > 0
                (3.0, 1.5, -1.2),  # e1 > 0; e2 < 0
            ],
            dtype=[("i_ixx", "<f8"), ("i_iyy", "<f8"), ("i_ixy", "<f8")],
        )

        cls.kwargs = {"band": "i"}

    def test_size(self):
        """Test CalcMomentSize functor"""
        traceSize = CalcMomentSize(sizeType="trace")(self.data, **self.kwargs)
        determinantSize = CalcMomentSize(sizeType="determinant")(self.data, **self.kwargs)

        for idx, row in enumerate(self.data):
            shape = Quadrupole(ixx=row["i_ixx"], iyy=row["i_iyy"], ixy=row["i_ixy"])
            self.assertFloatsAlmostEqual(traceSize[idx], shape.getTraceRadius(), rtol=1e-8)
            self.assertFloatsAlmostEqual(determinantSize[idx], shape.getDeterminantRadius(), rtol=1e-8)
            # Arithmetic mean >= Geometric mean implies that
            # trace radius is never smaller than determinant radius.
            self.assertGreaterEqual(traceSize[idx], determinantSize[idx])

    def test_complex_shear(self):
        """Test CalcE functor

        Test that our ellipticity calculation under the two conventions are
        accurate by comparing with GalSim routines.
        """
        shear = CalcE(ellipticityType="shear")(self.data, **self.kwargs)
        distortion = CalcE(ellipticityType="distortion")(self.data, **self.kwargs)
        size = CalcMomentSize(sizeType="determinant")(self.data, **self.kwargs)
        for idx, row in enumerate(self.data):
            galsim_shear = galsim.Shear(shear[idx])
            self.assertFloatsAlmostEqual(distortion[idx].real, galsim_shear.e1)
            self.assertFloatsAlmostEqual(distortion[idx].imag, galsim_shear.e2)
            # Check that the ellipiticity values correspond to the moments.
            A = galsim_shear.getMatrix() * size[idx]
            M = np.dot(A.transpose(), A)
            self.assertFloatsAlmostEqual(M[0, 0], row["i_ixx"], rtol=1e-8)
            self.assertFloatsAlmostEqual(M[1, 1], row["i_iyy"], rtol=1e-8)
            self.assertFloatsAlmostEqual(M[0, 1], row["i_ixy"], rtol=1e-8)

    def test_halve_angle(self):
        """Test ``halvePhaseAngle`` parameter in CalcE

        Test that setting ``halvePhaseAngle`` to True halves the phase angle
        while keeping the magnitude the same.
        """
        ellip = CalcE(ellipticityType="shear")(self.data, **self.kwargs)
        ellip_half = CalcE(ellipticityType="shear", halvePhaseAngle=True)(self.data, **self.kwargs)
        self.assertFloatsAlmostEqual(np.abs(ellip), np.abs(ellip_half))

        for idx, row in enumerate(self.data):
            galsim_shear = galsim.Shear(ellip[idx])
            galsim_shear_half = galsim.Shear(ellip_half[idx])
            self.assertFloatsAlmostEqual(np.abs(ellip_half[idx]), galsim_shear.g)
            self.assertFloatsAlmostEqual(
                galsim_shear.beta / galsim.radians, 2 * galsim_shear_half.beta / galsim.radians
            )

    @lsst.utils.tests.methodParameters(ellipticityType=("distortion", "shear"))
    def test_shear_components(self, ellipticityType):
        """Test CalcE1 and CalcE2 functors

        This test checks if CalcE1 and CalcE2 correspond to the real and
        imaginary components of CalcE.
        """
        ellip = CalcE(ellipticityType=ellipticityType)(self.data, **self.kwargs)
        e1 = CalcE1(ellipticityType=ellipticityType)(self.data, **self.kwargs)
        e2 = CalcE2(ellipticityType=ellipticityType)(self.data, **self.kwargs)

        self.assertFloatsAlmostEqual(np.real(ellip), e1)
        self.assertFloatsAlmostEqual(np.imag(ellip), e2)

    def test_e1_validation(self):
        """Test that CalcE1 throws an exception when misconfigured."""
        CalcE1(ellipticityType="distortion", colXy=None).validate()
        with self.assertRaises(FieldValidationError):
            CalcE1(ellipticityType="shear", colXy=None).validate()

    def test_ediff_validation(self):
        """Test that CalcEDiff takes ellipticities of same convention."""
        ellipA = CalcE(ellipticityType="shear")
        ellipB = CalcE(ellipticityType="distortion")
        with self.assertRaises(FieldValidationError):
            CalcEDiff(colA=ellipA, colB=ellipB).validate()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
