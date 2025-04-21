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


import os
import shutil
import tempfile
import unittest

import lsst.utils.tests
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from lsst.analysis.tools.actions.plot import CompletenessHist
from lsst.analysis.tools.actions.plot.plotUtils import get_and_remove_figure_text
from lsst.analysis.tools.math import divide, sqrt

# Set to True to debug plot test(s)
debug_plot = False
if not debug_plot:
    matplotlib.use("Agg")

ROOT = os.path.abspath(os.path.dirname(__file__))
filename_texts_ref = os.path.join(ROOT, "data", "test_completenessPlot_texts.txt")
path_lines_ref = os.path.join(ROOT, "data", "test_completenessPlot_lines")


class CompletenessPlotTestCase(lsst.utils.tests.TestCase):
    """CompletenessHist test case."""

    def setUp(self):
        self.testDir = tempfile.mkdtemp(dir=ROOT, prefix="test_output_completeness")

        plot = CompletenessHist()

        plot.action.action.selector_range_ref.vectorKey = "mag_ref"
        plot.action.action.selector_range_target.vectorKey = "mag_target"
        plot.action.bins.mag_low_max = 25500
        plot.action.bins.mag_interval = 100
        plot.action.bins.mag_width = 200

        band = "r"

        mag = 12.5 + 2.5 * np.log10(np.arange(10, 100000))
        zeropoint = mag[-1] + 1
        flux = 10 ** (-0.4 * (mag - zeropoint))
        rng = np.random.default_rng(0)
        # See usage below if changing this to allow nan/bad extendedness
        extendedness = 0.0 + (rng.uniform(size=len(mag)) < 0.99 * (mag - mag[0]) / (mag[-1] - mag[0]))
        flux_meas = flux + rng.normal(scale=sqrt(flux * (1 + extendedness)))
        flux_pos = flux_meas > 0
        flux_err = np.empty_like(flux_meas)
        flux_err[flux_pos] = sqrt(flux_meas[flux_pos])
        flux_err[~flux_pos] = 0
        detect_sn = divide(flux_meas, flux_err)
        good = detect_sn > 3
        match_distance = np.full_like(detect_sn, np.nan)
        rand_norm_abs = np.abs(rng.normal(size=match_distance.size))
        # make some unmatched as a function of S/N
        matched = rand_norm_abs < 0.2 * detect_sn
        match_distance[matched] = rand_norm_abs[matched]

        flux_meas[~good] = np.nan

        # estimated from a DC2 tract
        prob_ps = 0.5 * (1 + np.tanh((19.2 - mag) / 2.8))
        is_ps = rng.uniform(size=len(prob_ps)) < prob_ps
        is_extended = rng.uniform(size=len(prob_ps)) > prob_ps
        is_extended[~good] = np.nan

        # add some false detections
        false_positive = np.abs(rng.normal(size=match_distance.size)) > 0.5 * detect_sn
        flux_meas_extra = flux_meas[false_positive]
        nan_extra = np.full_like(flux_meas_extra, np.nan)
        action_hist = plot.action.action
        action_hist.key_mask_ref = "select_is_ref"
        action_hist.key_mask_target = "select_is_target"
        is_ps = np.concatenate((is_ps, nan_extra))
        extendedness = np.concatenate((extendedness, extendedness[false_positive]))
        is_extended = np.concatenate((is_extended, is_extended[false_positive]))

        data = {
            "mag_ref": np.concatenate((mag, nan_extra)),
            "mag_target": -2.5 * np.log10(np.concatenate((flux_meas, flux_meas_extra))) + zeropoint,
            "refcat_is_pointsource": is_ps,
            action_hist.key_match_distance: np.concatenate((match_distance, nan_extra)),
            action_hist.key_matched_class: is_ps == ~is_extended,
            action_hist.key_mask_ref: np.isfinite(is_ps),
            action_hist.key_mask_target: np.isfinite(extendedness),
        }

        results = plot.action(data, band="r")
        data.update(results)
        self.band = "r"
        self.data = data
        self.plot = plot
        self.plotInfo = {"plotName": "compurity test", "run": "", "tableName": "", "bands": band}

    def tearDown(self):
        if os.path.exists(self.testDir):
            shutil.rmtree(self.testDir, True)
        del self.data
        del self.plot
        del self.plotInfo
        del self.testDir

    def test_completenessPlotWithTwoHistsTask(self):
        plt.rcParams.update(plt.rcParamsDefault)
        result = self.plot(
            data=self.data,
            band=self.band,
            plotInfo=self.plotInfo,
        )
        self.assertTrue(isinstance(result, plt.Figure))

        to_compare = True
        if to_compare:
            # Set to true to save plots as PNGs
            # Use matplotlib.testing.compare.compare_images if needed
            save_images = False
            if save_images:
                result.savefig(os.path.join(ROOT, "data", "test_completenessPlot.png"))

            texts, lines = get_and_remove_figure_text(result)
            if save_images:
                result.savefig(os.path.join(ROOT, "data", "test_completenessPlot_unlabeled.png"))

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
                # nans are possible here (divide by zero)
                self.assertFloatsAlmostEqual(arr, line, atol=1e-3, rtol=1e-4, ignoreNaNs=True)

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

            self.assertEqual(texts_ref, texts_set)
        if debug_plot:
            plt.show()


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
