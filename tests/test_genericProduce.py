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

import lsst.utils.tests
import numpy as np
from lsst.analysis.tools.atools.genericProduce import ColumnMagnitudeScatterPlot


class ColumnMagnitudeScatterPlotTestCase(TestCase):
    """Test that analysis tools can be loaded from a pipeline"""

    def setUp(self) -> None:
        super().setUp()

        # Set up a quasi-plausible measurement catalog
        band = "i"
        mag = 12.5 + 2.5 * np.log10(np.arange(10, 100000))
        n_obj = len(mag)
        flux = 10 ** (-0.4 * (mag - (mag[-1] + 1)))
        rng = np.random.default_rng(0)
        extendedness = 0.0 + (rng.uniform(size=len(mag)) < 0.99 * (mag - mag[0]) / (mag[-1] - mag[0]))
        flux_meas = flux + rng.normal(scale=np.sqrt(flux * (1 + extendedness)))
        flux_err = np.sqrt(flux_meas * (1 + extendedness))
        good = (flux_meas / np.sqrt(flux * (1 + extendedness))) > 3
        extendedness = extendedness[good]
        flux_meas = flux_meas[good]
        flux_err = flux_err[good]
        key_flux = f"{band}_cModelFlux"

        data = {
            key_flux: flux_meas,
            f"{key_flux}Err": flux_err,
            "detect_isPrimary": np.ones(n_obj, dtype=bool),
            "refExtendedness": extendedness,
            "y": np.arange(n_obj),
        }

        self.action = ColumnMagnitudeScatterPlot(
            key_y="y",
            mag_x="cmodel_err",
        )
        self.action.finalize()
        self.band = band
        self.data = data  # Table(data)

    def testAction(self) -> None:
        results = self.action(self.data, band=self.band)
        assert tuple(results.keys()) == ("hist_NullPlot", "xy_Custom")


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
