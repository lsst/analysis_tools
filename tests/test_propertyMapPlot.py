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
import unittest

import healsparse as hsp
import lsst.utils.tests
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import skyproj
from lsst.analysis.tools.atools.healSparsePropertyMap import (
    PerTractPropertyMapTool,
    SurveyWidePropertyMapTool,
)
from lsst.analysis.tools.tasks.propertyMapAnalysis import (
    PerTractPropertyMapAnalysisConfig,
    PerTractPropertyMapAnalysisTask,
    SurveyWidePropertyMapAnalysisConfig,
    SurveyWidePropertyMapAnalysisTask,
)
from lsst.daf.butler import Butler, DataCoordinate, DatasetType, DeferredDatasetHandle
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.skymap.discreteSkyMap import DiscreteSkyMap
from mpl_toolkits import axisartist

# No display needed.
matplotlib.use("Agg")

# Direcory where this file is located.
ROOT = os.path.abspath(os.path.dirname(__file__))


class PerTractPropertyMapAnalysisTaskTestCase(lsst.utils.tests.TestCase):
    """PerTractPropertyMapAnalysisTask test case.

    Notes
    -----
    While definitive tests are conducted in `ci_hsc` and `ci_imsim` using real
    and simulated datasets to ensure thorough coverage, this test case is
    designed to catch foundational issues like syntax errors or logical
    inconsistencies in the way the plots are generated.
    """

    def setUp(self):
        # Create a temporary directory to test in.
        self.testDir = makeTestTempDir(ROOT)

        # Create a butler in the test directory.
        Butler.makeRepo(self.testDir)
        butler = Butler(self.testDir, run="testrun")

        # Make a dummy dataId.
        dataId = {"band": "i", "skymap": "hsc_rings_v1", "tract": 1915}
        dataId = DataCoordinate.standardize(dataId, universe=butler.dimensions)

        # Configure the maps to be plotted.
        config = PerTractPropertyMapAnalysisConfig()

        # Set configurations sent to skyproj.
        config.projectionKwargs = {"celestial": True, "gridlines": True}
        config.colorbarKwargs = {"location": "top", "cmap": "viridis"}

        # The entries in the 'atools' namespace must exactly match the dataset
        # types.
        config.atools.deepCoadd_exposure_time_map_sum = PerTractPropertyMapTool()
        config.atools.deepCoadd_exposure_time_map_sum.nBinsHist = 24
        config.atools.deepCoadd_psf_maglim_map_weighted_mean = PerTractPropertyMapTool()
        config.atools.goodSeeingCoadd_dcr_dra_map_weighted_mean = PerTractPropertyMapTool()

        # Generate a list of dataset type names.
        names = [name for name in config.atools.fieldNames]

        # Mock up corresponding HealSparseMaps and register them with the
        # butler.
        inputs = {}
        for name, value in zip(names, np.linspace(1, 10, len(names))):
            hspMap = hsp.HealSparseMap.make_empty(nside_coverage=32, nside_sparse=4096, dtype=np.float32)
            hspMap[0:10000] = value
            hspMap[100000:110000] = value + 1
            hspMap[500000:510000] = value + 2
            datasetType = DatasetType(name, [], "HealSparseMap", universe=butler.dimensions)
            butler.registry.registerDatasetType(datasetType)
            dataRef = butler.put(hspMap, datasetType)
            # Keys in inputs are designed to reflect the task's connection
            # names.
            inputs[name] = DeferredDatasetHandle(butler=butler, ref=dataRef, parameters=None)

        # Mock up the skymap and tractInfo.
        skyMapConfig = DiscreteSkyMap.ConfigClass()
        coords = [  # From the PS1 Medium-Deep fields.
            (10.6750, 41.2667),  # M31
            (36.2074, -04.5833),  # XMM-LSS
        ]
        skyMapConfig.raList = [c[0] for c in coords]
        skyMapConfig.decList = [c[1] for c in coords]
        skyMapConfig.radiusList = [2] * len(coords)
        skyMapConfig.validate()
        skymap = DiscreteSkyMap(config=skyMapConfig)
        self.tractInfo = skymap.generateTract(0)

        # Initialize the task and set class attributes for subsequent use.
        task = PerTractPropertyMapAnalysisTask()
        self.config = config
        self.plotInfo = task.parsePlotInfo(inputs, dataId, list(inputs.keys()))
        self.data = inputs

        for atool in self.config.atools:
            atool.finalize()

    def tearDown(self):
        del self.data
        del self.config
        del self.tractInfo
        del self.plotInfo
        removeTestTempDir(self.testDir)
        del self.testDir

    def test_PerTractPropertyMapAnalysisTask(self):
        # Previously computed reference RGB fractions for the plots.
        expectedRGBFractions = [
            (0.3195957485702613, 0.3214485564746734, 0.3205870925245099),
            (0.3196810334967318, 0.3215274948937909, 0.32067869434232027),
            (0.31861285794526134, 0.3206199550653596, 0.3196505213439543),
            (0.3185880596405229, 0.3206015032679739, 0.319619406147876),
            (0.31843414266748354, 0.32044758629493475, 0.3194654891748366),
        ]

        plt.rcParams.update(plt.rcParamsDefault)
        for name, atool, expectedRGBFraction in zip(
            self.config.atools.fieldNames, self.config.atools, expectedRGBFractions
        ):
            # Run the task via butler using the tool.
            result = atool(
                data=self.data, tractInfo=self.tractInfo, plotConfig=self.config, plotInfo=self.plotInfo
            )
            key = atool.process.buildActions.data.mapKey + "_PerTractPropertyMapPlot"

            # Check that the key is in the result.
            self.assertIn(key, result)

            # Check that the output is a matplotlib figure.
            fig = result[key]
            self.assertTrue(isinstance(fig, plt.Figure), msg=f"Figure {key} is not a matplotlib figure.")

            # Assert the number of axes in the figure. At least not empty.
            self.assertEqual(len(fig.axes), 10)

            # Verify some configuration parameters.
            binsCount = atool.nBinsHist
            if "exposure" in name.lower():
                self.assertEqual(binsCount, 24)
            else:
                self.assertEqual(binsCount, 100)
            zoomFactors = self.config.zoomFactors
            self.assertEqual(len(zoomFactors), 2)

            # Extract the property name from the dataset type name and infer
            # the x-axis label of the histogram.
            propertyName = name.split("_map_")[0].split("Coadd_")[1].replace("_", " ")
            xlabel = (
                propertyName.title()
                .replace("Psf", "PSF")
                .replace("Dcr", "DCR")
                .replace("Dra", r"$\Delta$RA")
                .replace("Ddec", r"$\Delta$Dec")
                .replace("E1", "e1")
                .replace("E2", "e2")
            )

            # Validate the structure of the figure.
            self._validateFigureStructure(fig, atool, binsCount, xlabel, zoomFactors)

            # Validate the RGB fractions of the figure. The tolerance is set
            # empirically.
            self._validateRGBFractions(fig, expectedRGBFraction, rtol=5e-3)

    @staticmethod
    def _isHistogramAxis(ax, binsCount, legendLabels, errors):
        """Checks if a given axis is a histogram axis based on specified
        parameters.

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            The axis to be checked.
        binsCount : `int`
            The expected number of bins in the histogram.
        legendLabels : `List` [`str`]
            The expected labels in the histogram legend.
        errors : `List` [`str`]
            A list to append any errors found during the checks.

        Returns
        -------
        None
            Errors are appended to the provided `errors` list.
        """

        # Count rectangle and polygon patches.
        nRectanglePatches = sum(1 for patch in ax.patches if isinstance(patch, matplotlib.patches.Rectangle))
        nPolygonPatches = sum(1 for patch in ax.patches if isinstance(patch, matplotlib.patches.Polygon))

        # Check for the number of rectangle patches for the filled histogram.
        if nRectanglePatches != binsCount:
            errors.append(
                f"Expected {binsCount} rectangle patches in histogram, but found {nRectanglePatches}."
            )

        # Check for the number of polygon patches, i.e. the step histograms.
        if nPolygonPatches != 2:
            errors.append(f"Expected 2 polygon patches, but found {nPolygonPatches}.")

        # Check for `fill_between` regions, represented by `PolyCollection`
        # objects.
        if len(ax.collections) != 2:
            errors.append(f"Expected 2 `fill_between` regions but found {len(ax.collections)}.")

        # Verify legend labels.
        legend = ax.get_legend()
        if not legend:
            errors.append("Legend is missing in the histogram.")
        else:
            labels = [text.get_text() for text in legend.get_texts()]
            if set(labels) != set(legendLabels):
                errors.append(f"Expected legend labels {legendLabels} but found {labels}.")

    def _validateFigureStructure(self, fig, atool, binsCount, xlabel, zoomFactors):
        """Validates the structure of a given matplotlib figure generated by
        the tool that is being tested.

        Parameters
        ----------
        fig : `~matplotlib.figure.Figure`
            The figure to be validated.
        atool :
            `~lsst.analysis.tools.atools.healSparsePropertyMap.
            PerTractPropertyMapTool`
            The tool that generated the figure.
        binsCount : `int`
            The expected number of bins in the histogram.
        xlabel : `str`
            The expected x-axis label of the histogram.
        zoomFactors : `List` [`float`]
            A list of zoom factors used for the zoomed-in plots.

        Raises
        ------
        AssertionError
            If any of the criteria for figure structure is not met. The error
            message will list all criteria that were not satisfied.
        """
        errors = []
        axes = fig.get_axes()

        # Check the total number of each axis type.
        totalSkyAxes = sum(isinstance(ax, skyproj.skyaxes.SkyAxes) for ax in axes)
        totalAxisArtistAxes = sum(isinstance(ax, axisartist.axislines.Axes) for ax in axes)
        totalColorbarAxes = sum(isinstance(ax, plt.Axes) and ax.get_label() == "<colorbar>" for ax in axes)

        if totalSkyAxes != 3:
            errors.append(f"Expected 3 SkyAxes but got {totalSkyAxes}.")
        if totalAxisArtistAxes != 3:
            errors.append(f"Expected 3 AxisArtist Axes but got {totalAxisArtistAxes}.")
        if totalColorbarAxes != 3:
            errors.append(f"Expected 3 colorbar Axes but got {totalColorbarAxes}.")

        # Check histogram axis.
        self._isHistogramAxis(
            axes[0],
            binsCount,
            ["Full Tract"]
            + [f"{atool.produce.plot.prettyPrintFloat(factor)}x Zoom" for factor in zoomFactors],
            errors,
        )

        # Verify x and y labels for histogram.
        if axes[0].get_xlabel() != xlabel:
            errors.append(f"Expected x-label '{xlabel}' for histogram but found '{axes[0].get_xlabel()}'.")
        if axes[0].get_ylabel() != "Normalized Count":
            errors.append(
                f"Expected y-label 'Normalized Count' for histogram but found '{axes[0].get_ylabel()}'."
            )

        self.assertTrue(len(errors) == 0, msg="\n" + "\n".join(errors))

    def _validateRGBFractions(self, fig, RGBFraction, rtol=1e-7):
        """Checks if a matplotlib figure has specified fractions of R, G, and B
        colors.

        Parameters
        ----------
        fig : `~matplotlib.figure.Figure`
            The figure to check.
        RGBFraction : `tuple`
            Tuple containing the desired fractions for red, green, and blue in
            the image, respectively.
        rtol : `float`, optional
            The relative tolerance allowed for the fractions. Default is 1e-7.

        Raises
        ------
        AssertionError
            If the actual fractions of the RGB colors in the image do not match
            the expected fractions within the given tolerance.
        """

        # Unpack the desired fractions.
        rFraction, gFraction, bFraction = RGBFraction

        # Draw the figure so the renderer can grab the pixel buffer.
        fig.canvas.draw()

        # Convert figure to data array.
        data = np.array(fig.canvas.renderer.buffer_rgba())[:, :, :3] / 255.0

        # Calculate fractions.
        rActualFraction = np.sum(data[:, :, 0]) / data.size
        gActualFraction = np.sum(data[:, :, 1]) / data.size
        bActualFraction = np.sum(data[:, :, 2]) / data.size

        # Check if the actual fractions meet the expected fractions within the
        # given tolerance.
        errors = []
        if not np.abs(rActualFraction - rFraction) <= rtol:
            errors.append(
                f"Calculated red fraction {rActualFraction} does not match {rFraction} within rtol {rtol}."
            )

        if not np.abs(gActualFraction - gFraction) <= rtol:
            errors.append(
                f"Calculated green fraction {gActualFraction} does not match {gFraction} within rtol {rtol}."
            )

        if not np.abs(bActualFraction - bFraction) <= rtol:
            errors.append(
                f"Calculated blue fraction {bActualFraction} does not match {bFraction} within rtol {rtol}."
            )

        self.assertTrue(len(errors) == 0, msg="\n" + "\n".join(errors))


class SurveyWidePropertyMapAnalysisTaskTestCase(lsst.utils.tests.TestCase):
    """PerTractPropertyMapAnalysisTask test case.

    Notes
    -----
    This is a basic functionality test to verify the internal workings of the
    task.
    """

    def setUp(self):
        # Create a temporary directory to test in.
        self.testDir = makeTestTempDir(ROOT)

        # Create a butler in the test directory.
        Butler.makeRepo(self.testDir)
        butler = Butler(self.testDir, run="testrun")

        # Make a dummy dataId.
        dataId = {"band": "i", "skymap": "hsc_rings_v1"}
        dataId = DataCoordinate.standardize(dataId, universe=butler.dimensions)

        # Configure the maps to be plotted.
        config = SurveyWidePropertyMapAnalysisConfig()

        # Set configurations sent to skyproj.
        config.autozoom = True
        config.projection = "Mollweide"
        config.projectionKwargs = {"celestial": True, "gridlines": True, "lon_0": 0}
        config.colorbarKwargs = {"location": "top", "cmap": "viridis"}

        # The entries in the 'atools' namespace must exactly match the dataset
        # type.
        config.atools.deepCoadd_exposure_time_consolidated_map_sum = SurveyWidePropertyMapTool()
        config.atools.deepCoadd_psf_maglim_consolidated_map_weighted_mean = SurveyWidePropertyMapTool()
        config.atools.goodSeeingCoadd_dcr_dra_consolidated_map_weighted_mean = SurveyWidePropertyMapTool()

        # Generate a list of dataset type names.
        names = [name for name in config.atools.fieldNames]

        # Mock up corresponding HealSparseMaps and register them with the
        # butler.
        inputs = {}
        for name, value in zip(names, np.linspace(1, 10, len(names))):
            hspMap = hsp.HealSparseMap.make_empty(nside_coverage=32, nside_sparse=4096, dtype=np.float32)
            hspMap[0:10000] = value
            hspMap[100000:110000] = value + 1
            hspMap[500000:510000] = value + 2
            datasetType = DatasetType(name, [], "HealSparseMap", universe=butler.dimensions)
            butler.registry.registerDatasetType(datasetType)
            dataRef = butler.put(hspMap, datasetType)
            # Keys in inputs are designed to reflect the dataset type names.
            inputs[name] = DeferredDatasetHandle(butler=butler, ref=dataRef, parameters=None)

        # Initialize the task and set class attributes for subsequent use.
        task = SurveyWidePropertyMapAnalysisTask()
        self.config = config
        self.plotInfo = task.parsePlotInfo(inputs, dataId, list(inputs.keys()))
        self.data = inputs

        for atool in self.config.atools:
            atool.finalize()

    def tearDown(self):
        del self.data
        del self.config
        del self.plotInfo
        removeTestTempDir(self.testDir)
        del self.testDir

    def test_SurveyWidePropertyMapAnalysisTask(self):
        plt.rcParams.update(plt.rcParamsDefault)
        for atool in self.config.atools:
            # Run the task via butler using the tool.
            result = atool(data=self.data, plotConfig=self.config, plotInfo=self.plotInfo)
            key = atool.process.buildActions.data.mapKey + "_SurveyWidePropertyMapPlot"

            # Check that the key is in the result.
            self.assertIn(key, result)

            # Check that the output is a matplotlib figure.
            fig = result[key]
            self.assertTrue(isinstance(fig, plt.Figure), msg=f"Figure {key} is not a matplotlib figure.")

            # Assert the number of axes in the figure. At least not empty.
            self.assertEqual(len(fig.axes), 3)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
