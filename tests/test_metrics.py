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
        kwargs = {key: "" for key in ("unit", "name_prefix")}
        for kwarg in kwargs:
            kwargs_init = {kwarg: ""}
            tester = MatchedRefCoaddMetric(**kwargs_init)
            # Validation will fail because configureMetrics hasn't been called
            with self.assertRaises(ValueError):
                tester.validate()
            # Failing to find any of the required keys
            with self.assertRaises(KeyError):
                tester({})
        tester = MatchedRefCoaddMetric(**kwargs)
        with self.assertRaises(ValueError):
            tester.validate()

        tester.finalize()
        # This works, although it leaves none keys in the schema
        # Should maybe be caught in populatePrepFromProcess
        inputs = list(tester.getInputSchema())
        n_input = len(inputs)
        self.assertGreater(n_input, 0)
        self.assertGreater(len(list(tester.configureMetrics())), 0)

        tester.finalize()
        # No more None key
        self.assertEquals(len(list(tester.getInputSchema())), n_input - 1)

    def testMatchedRefCoaddDiffMagMetric(self):
        self._testMatchedRefCoaddMetricDerived(
            MatchedRefCoaddDiffMagMetric,
            fluxes={"cmodel": MagnitudeTool.fluxes_default.cmodel_err},
            mag_y="cmodel",
            name_prefix="",
            unit="",
        )

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
