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
from lsst.analysis.tools.actions.vector import ConstantValue, MultiplyVector
from lsst.analysis.tools.atools import (
    MagnitudeTool,
    MatchedRefCoaddCompurityTool,
    MatchedRefCoaddDiffColorTool,
    MatchedRefCoaddDiffDistanceTool,
    MatchedRefCoaddDiffMagTool,
    MatchedRefCoaddDiffPositionTool,
    MatchedRefCoaddDiffTool,
)


class MatchedRefCoaddDiffToolMinimal(MatchedRefCoaddDiffTool):
    """A bare-minimum implementation of a diff tool."""

    def get_key_flux_y(self) -> str:
        return self.mag_x

    def finalize(self):
        if not hasattr(self.process.buildActions, "diff"):
            super().finalize()
            self._set_actions()

            self.process.buildActions.diff = MultiplyVector(
                actionA=getattr(self.process.buildActions, f"mag_{self.mag_x}"),
                actionB=ConstantValue(value=0.001),
            )


class TestDiffMatched(TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.band_default = "analysisTools"

    def _testMatchedRefCoaddMetricDerived(
        self,
        type_metric: type[MatchedRefCoaddCompurityTool | MatchedRefCoaddDiffColorTool],
        suffixes_configure: list[str] | None = None,
        **kwargs,
    ):
        if suffixes_configure is None:
            suffixes_configure = [""]
        plotInfo = {key: "" for key in ("plotName", "run", "tableName")}
        plotInfo["bands"] = []

        tester = type_metric(**kwargs)
        has_compute_chi = hasattr(tester, "compute_chi")
        has_configureMetrics = hasattr(tester, "configureMetrics")

        for compute_chi in (False, True) if has_compute_chi else (None,):
            if has_compute_chi:
                tester.compute_chi = compute_chi
            # tester.getInputSchema won't work properly before finalizing
            tester.finalize()
            keys = list({k[0]: None for k in tester.getInputSchema()})
            self.assertGreater(len(keys), 0)
            if has_configureMetrics:
                for suffix in suffixes_configure:
                    self.assertGreater(len(list(tester.configureMetrics(attr_suffix=suffix))), 0)
            data = {}
            n_data = 10
            for key in keys:
                if key.endswith("is_pointsource"):
                    values = (np.arange(0, n_data) % 2) == 1
                else:
                    values = np.linspace(0.1, 15.0, n_data)
                data[key.format(band=self.band_default)] = values

            output = tester(data, skymap=None, plotInfo=plotInfo)
            self.assertGreater(len(output), 0)
            if has_compute_chi:
                tester = type_metric(**kwargs)

    def testMatchedRefCoaddMetric(self):
        kwargs = {key: "" for key in ("unit", "name_prefix")}
        # The metric can now be set up with default kwargs
        # Pass one at a time to test
        for kwarg in kwargs:
            kwargs_init = {kwarg: ""}
            tester = MatchedRefCoaddDiffToolMinimal(**kwargs_init)
            tester.validate()
            with self.assertRaises(KeyError):
                tester.configureMetrics()
            tester.finalize()
            tester.configureMetrics()
            # Failing to find any of the required keys
            with self.assertRaises(KeyError):
                tester({})
        tester = MatchedRefCoaddDiffToolMinimal(**kwargs)
        tester.validate()
        with self.assertRaises(KeyError):
            tester.configureMetrics()
        tester.finalize()
        inputs = list(tester.getInputSchema())
        n_input = len(inputs)
        self.assertGreater(n_input, 0)
        self.assertGreater(len(list(tester.configureMetrics())), 0)

        self.assertEqual(len(inputs), n_input)
        data = {key.format(band="analysisTools"): np.array([0.0]) for key, *_ in inputs}
        # There's no metric or plot so it just returns an empty dict
        self.assertEqual(len(tester(data)), 0)

    def testMatchedRefCoaddCompurity(self):
        self._testMatchedRefCoaddMetricDerived(
            MatchedRefCoaddCompurityTool,
        )

    def testMatchedRefCoaddDiffColor(self):
        self._testMatchedRefCoaddMetricDerived(
            MatchedRefCoaddDiffColorTool,
            suffixes_configure=["_0", "_1"],
            fluxes={"cmodel": MagnitudeTool.fluxes_default.cmodel_err},
            mag_y1="cmodel_err",
            mag_y2="cmodel_err",
            bands={
                "g": "i",
                "u": "z",
            },
            name_prefix="",
            unit="",
        )

    def testMatchedRefCoaddDiffDistance(self):
        self._testMatchedRefCoaddMetricDerived(
            MatchedRefCoaddDiffDistanceTool,
        )

    def testMatchedRefCoaddDiffMag(self):
        self._testMatchedRefCoaddMetricDerived(
            MatchedRefCoaddDiffMagTool,
            fluxes={"cmodel": MagnitudeTool.fluxes_default.cmodel_err},
            mag_y="cmodel",
            name_prefix="",
            unit="",
        )

    def testMatchedRefCoaddDiffPosition(self):
        for variable in ("x", "y"):
            self._testMatchedRefCoaddMetricDerived(
                MatchedRefCoaddDiffPositionTool,
                coord_meas=variable,
                coord_ref=variable,
            )


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
