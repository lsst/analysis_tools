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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from unittest import TestCase, main

import lsst.pex.config as pexConfig
import lsst.utils.tests

# This is needed to work around a bug with pyproj when it is imported
# inside a function, in this case through a pipeline
import skyproj  # noqa: F401
from lsst.analysis.tools.interfaces import AnalysisTool


class NamedTool(AnalysisTool):
    name = pexConfig.Field[str](doc="Name", default="")


class A(NamedTool):
    pass


class B(NamedTool):
    def finalize(self):
        super().finalize()
        self.name += "B"


class C(B):
    # Override the default finalize to show that it does get called just once
    # at the end. This should not be implemented normally.
    def _baseFinalize(self) -> None:
        self.name += "_final"


class D(A, C):
    def finalize(self):
        super().finalize()
        self.name += "D"


class E(C, A):
    def finalize(self):
        super().finalize()
        self.name += "E"


class FinalizeTestCase(TestCase):
    """Test that finalize method properly calls through the MRO, and that
    it appropriately calls the base finalize once, at the end.
    """

    def setUp(self) -> None:
        super().setUp()
        self.a = A()
        self.c = C()
        self.d = D()
        self.e = E()

    def testFinalize(self):
        self.a.finalize()
        self.assertEqual(self.a.name, "")
        self.c.finalize()
        self.assertEqual(self.c.name, "B_final")
        self.d.finalize()
        self.assertEqual(self.d.name, "BD_final")
        self.e.finalize()
        self.assertEqual(self.e.name, "BE_final")


class TestMetricNameValidation(TestCase):
    """Test that disallowed characters in measurement names are filtered."""

    def testMinusNameFails(self):
        tool = AnalysisTool()
        with self.assertRaises(pexConfig.FieldValidationError):
            tool.produce.metric.newNames = {"something": "not-allowed"}


class LoadFromPipelineTestCase(TestCase):
    """Test that analysis tools can be loaded from a pipeline"""

    def setUp(self) -> None:
        super().setUp()
        self.pipeline = "$ANALYSIS_TOOLS_DIR/pipelines/coaddQualityExtended.yaml"

    def testLoadRelative(self) -> None:
        # test that relative names can be loaded
        tool = AnalysisTool.fromPipeline(self.pipeline, "shapeSizeDetRadiusVsCmodelMag")
        self.assertIsNotNone(tool)
        self.assertEqual(tool.mag_x, "cmodel_err")  # type: ignore
        self.assertEqual(tool.size_y, "shape_slot")  # type: ignore

    def testLoadAbs(self) -> None:
        ...
        # test that relative names can be loaded
        tool = AnalysisTool.fromPipeline(self.pipeline, "atools.shapeSizeDetRadiusVsCmodelMag", fullpath=True)
        self.assertIsNotNone(tool)
        self.assertEqual(tool.mag_x, "cmodel_err")  # type: ignore
        self.assertEqual(tool.size_y, "shape_slot")  # type: ignore


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
