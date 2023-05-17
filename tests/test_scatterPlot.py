# This file is part of analysis_drp.
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


import os
import shutil
import tempfile
import unittest

import lsst.utils.tests
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lsst.analysis.tools.actions.plot.plotUtils import get_and_remove_figure_text
from lsst.analysis.tools.actions.plot.scatterplotWithTwoHists import (
    ScatterPlotStatsAction,
    ScatterPlotWithTwoHists,
)
from lsst.analysis.tools.actions.vector.mathActions import ConstantValue, DivideVector, SubtractVector
from lsst.analysis.tools.actions.vector.selectors import (
    GalaxySelector,
    SnSelector,
    StarSelector,
    VectorSelector,
)
from lsst.analysis.tools.actions.vector.vectorActions import ConvertFluxToMag, DownselectVector, LoadVector
from lsst.analysis.tools.interfaces import AnalysisTool

matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.dirname(__file__))
filename_texts_ref = os.path.join(ROOT, "data", "test_scatterPlot_texts.txt")
path_lines_ref = os.path.join(ROOT, "data", "test_scatterPlot_lines")


class ScatterPlotWithTwoHistsTaskTestCase(lsst.utils.tests.TestCase):
    """ScatterPlotWithTwoHistsTask test case."""

    def setUp(self):
        self.testDir = tempfile.mkdtemp(dir=ROOT, prefix="test_output")

        # Set up a quasi-plausible measurement catalog
        mag = 12.5 + 2.5 * np.log10(np.arange(10, 100000))
        flux = 10 ** (-0.4 * (mag - (mag[-1] + 1)))
        rng = np.random.default_rng(0)
        extendedness = 0.0 + (rng.uniform(size=len(mag)) < 0.99 * (mag - mag[0]) / (mag[-1] - mag[0]))
        flux_meas = flux + rng.normal(scale=np.sqrt(flux * (1 + extendedness)))
        flux_err = np.sqrt(flux_meas * (1 + extendedness))
        good = (flux_meas / np.sqrt(flux * (1 + extendedness))) > 3
        extendedness = extendedness[good]
        flux = flux[good]
        flux_meas = flux_meas[good]
        flux_err = flux_err[good]

        # Configure the plot to show observed vs true mags
        action = ScatterPlotWithTwoHists(
            xAxisLabel="mag",
            yAxisLabel="mag meas - ref",
            magLabel="mag",
            plotTypes=[
                "galaxies",
                "stars",
            ],
            xLims=(20, 30),
            yLims=(-1000, 1000),
        )
        plot = AnalysisTool()
        plot.produce.plot = action

        # Load the relevant columns
        key_flux = "meas_Flux"
        plot.process.buildActions.fluxes_meas = LoadVector(vectorKey=key_flux)
        plot.process.buildActions.fluxes_err = LoadVector(vectorKey=f"{key_flux}Err")
        plot.process.buildActions.fluxes_ref = LoadVector(vectorKey="ref_Flux")
        plot.process.buildActions.mags_ref = ConvertFluxToMag(
            vectorKey=plot.process.buildActions.fluxes_ref.vectorKey
        )

        # Compute the y-axis quantity
        plot.process.buildActions.diff = SubtractVector(
            actionA=ConvertFluxToMag(
                vectorKey=plot.process.buildActions.fluxes_meas.vectorKey, returnMillimags=True
            ),
            actionB=DivideVector(
                actionA=plot.process.buildActions.mags_ref,
                actionB=ConstantValue(value=1e-3),
            ),
        )

        # Filter stars/galaxies, storing quantities separately
        plot.process.buildActions.galaxySelector = GalaxySelector(vectorKey="refExtendedness")
        plot.process.buildActions.starSelector = StarSelector(vectorKey="refExtendedness")
        for singular, plural in (("galaxy", "Galaxies"), ("star", "Stars")):
            setattr(
                plot.process.filterActions,
                f"x{plural}",
                DownselectVector(
                    vectorKey="mags_ref", selector=VectorSelector(vectorKey=f"{singular}Selector")
                ),
            )
            setattr(
                plot.process.filterActions,
                f"y{plural}",
                DownselectVector(vectorKey="diff", selector=VectorSelector(vectorKey=f"{singular}Selector")),
            )
            setattr(
                plot.process.filterActions,
                f"flux{plural}",
                DownselectVector(
                    vectorKey="fluxes_meas", selector=VectorSelector(vectorKey=f"{singular}Selector")
                ),
            )
            setattr(
                plot.process.filterActions,
                f"flux{plural}Err",
                DownselectVector(
                    vectorKey="fluxes_err", selector=VectorSelector(vectorKey=f"{singular}Selector")
                ),
            )

            # Compute low/high SN summary stats
            statAction = ScatterPlotStatsAction(
                vectorKey=f"y{plural}",
                fluxType=f"flux{plural}",
                highSNSelector=SnSelector(fluxType=f"flux{plural}", threshold=50),
                lowSNSelector=SnSelector(fluxType=f"flux{plural}", threshold=20),
            )
            setattr(plot.process.calculateActions, plural.lower(), statAction)

        data = {
            "ref_Flux": flux,
            key_flux: flux_meas,
            f"{key_flux}Err": flux_err,
            "refExtendedness": extendedness,
        }

        self.data = pd.DataFrame(data)
        self.plot = plot
        self.plot.finalize()
        plotInfo = {key: "test" for key in ("plotName", "run", "tableName")}
        plotInfo["bands"] = []
        self.plotInfo = plotInfo

    def tearDown(self):
        if os.path.exists(self.testDir):
            shutil.rmtree(self.testDir, True)
        del self.data
        del self.plot
        del self.plotInfo
        del self.testDir

    def test_ScatterPlotWithTwoHistsTask(self):
        plt.rcParams.update(plt.rcParamsDefault)
        result = self.plot(
            data=self.data,
            skymap=None,
            plotInfo=self.plotInfo,
        )
        # unpack the result from the dictionary
        result = result[type(self.plot.produce.plot).__name__]
        self.assertTrue(isinstance(result, plt.Figure))

        # Set to true to save plots as PNGs
        # Use matplotlib.testing.compare.compare_images if needed
        save_images = False
        if save_images:
            result.savefig(os.path.join(ROOT, "data", "test_scatterPlot.png"))

        texts, lines = get_and_remove_figure_text(result)
        if save_images:
            result.savefig(os.path.join(ROOT, "data", "test_scatterPlot_unlabeled.png"))

        # Set to true to re-generate reference data
        resave = False

        # Compare line values
        for idx, line in enumerate(lines):
            filename = os.path.join(path_lines_ref, f"line_{idx}.txt")
            if resave:
                np.savetxt(filename, line)
            arr = np.loadtxt(filename)
            # Differences of order 1e-12 possible between MacOS and Linux
            # Plots are generally not expected to be that precise
            # Differences to 1e-3 should not be visible with this test data
            self.assertFloatsAlmostEqual(arr, line, atol=1e-3, rtol=1e-4)

        # Ensure that newlines within labels are replaced by a sentinel
        newline = "\n"
        newline_replace = "[newline]"
        # Compare text labels
        if resave:
            with open(filename_texts_ref, "w") as f:
                f.writelines(f"{text.strip().replace(newline, newline_replace)}\n" for text in texts)

        with open(filename_texts_ref, "r") as f:
            texts_ref = set(x.strip() for x in f.readlines())
        texts_set = set(x.strip().replace(newline, newline_replace) for x in texts)

        self.assertTrue(texts_ref.issuperset(texts_set))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
