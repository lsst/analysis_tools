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
    pass


class D(A, C):
    def finalize(self):
        super().finalize()
        self.name += "D"


class E(C, A):
    def finalize(self):
        super().finalize()
        self.name += "E"


class FinalizeTestCase(TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.a = A()
        self.c = C()
        self.d = D()
        self.e = E()

    def testFinalize(self):
        self.a.finalize()
        self.assertEquals(self.a.name, "")
        self.c.finalize()
        self.assertEquals(self.c.name, "B")
        self.d.finalize()
        self.assertEquals(self.d.name, "BD")
        self.e.finalize()
        self.assertEquals(self.e.name, "BE")


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
