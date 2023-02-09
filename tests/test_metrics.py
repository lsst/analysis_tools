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
from __future__ import annotations

from unittest import TestCase, main

import lsst.utils.tests
import numpy as np
from lsst.analysis.tools.atools import (
    MatchedRefCoaddCModelFluxMetric,
    MatchedRefCoaddMetric,
    MatchedRefCoaddPositionMetric,
)
from lsst.analysis.tools.contexts import MatchedRefChiContext, MatchedRefDiffContext


class TestDiffMatched(TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.band_default = "analysisTools"
        self.contexts = (MatchedRefChiContext, MatchedRefDiffContext)

    def _testMatchedRefCoaddMetricDerived(self, type_metric: type[MatchedRefCoaddMetric], **kwargs):
        for context in self.contexts:
            tester = type_metric(**kwargs)
            tester.applyContext(context)
            tester.finalize()

            keys = list(k[0] for k in tester.getInputSchema())
            self.assertGreater(len(keys), 0)
            self.assertGreater(len(list(tester.configureMetrics())), 0)
            data = {key.format(band=self.band_default): np.arange(5) for key in keys}
            self.assertGreater(len(tester(data)), 0)

    def testMatchedRefCoaddMetric(self):
        tester = MatchedRefCoaddMetric(unit="")
        with self.assertRaises(ValueError):
            tester({})
        tester = MatchedRefCoaddMetric(name_prefix="")
        with self.assertRaises(ValueError):
            tester({})
        tester = MatchedRefCoaddMetric(unit="", name_prefix="")
        tester.finalize()
        self.assertGreater(len(list(tester.getInputSchema())), 0)
        self.assertGreater(len(list(tester.configureMetrics())), 0)

    def testMatchedRefCoaddCModelFluxMetric(self):
        self._testMatchedRefCoaddMetricDerived(MatchedRefCoaddCModelFluxMetric)

    def testMatchedRefCoaddPositionMetric(self):
        for variable in ("x", "y"):
            self._testMatchedRefCoaddMetricDerived(MatchedRefCoaddPositionMetric, variable=variable)


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
