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
from __future__ import annotations

__all__ = ("WholeSkyPlotTool",)

import re
from string import Formatter

from lsst.pex.config import Field, ListField

from ..actions.plot.wholeSkyPlot import WholeSkyPlot
from ..interfaces import AnalysisTool


class WholeSkyPlotTool(AnalysisTool):
    """Plot metrics across all tracts on the sky."""

    parameterizedBand = False
    propagateData = True

    bands = ListField[str](
        doc="Photometric bands to use for plots.",
        default=["u", "g", "r", "i", "z", "y"],
        listCheck=lambda x: len(set(x)) == len(x),
    )

    plotKeys = ListField[str](
        doc="Metrics to plot that are not band-specific.",
        default=["yPerpPSF_yPerp_psfFlux_median"],
        listCheck=lambda x: len(set(x)) == len(x),
    )
    keysWithBand = ListField[str](
        doc="Metrics to plot that are band-specific.",
        default=["e1Diff_{band}_highSNStars_median"],
        itemCheck=lambda x: "{band}" in x,
        listCheck=lambda x: len(set(x)) == len(x),
    )
    xAxisLabel = Field[str](doc="Label for the x axis.", default="RA (degrees)")
    yAxisLabel = Field[str](doc="Label for the y axis.", default="Dec (degrees)")
    autoAxesLimits = Field[bool](doc="Find axes limits automatically.", default=True)
    xLimits = ListField[float](doc="Plotting limits for the x axis.", default=[-5.0, 365.0])
    yLimits = ListField[float](doc="Plotting limits for the y axis.", default=[-10.0, 60.0])
    figureSize = ListField[float](doc="Size of the figure.", default=[9.0, 3.5])
    colorBarRange = Field[float](
        doc="The multiplier for the color bar range. The max/min range values are: median +/- N * sigmaMad"
        ", where N is this config value.",
        default=3.0,
    )
    sequentialMetrics = ListField[str](
        doc="Partial names of metrics with sequential values. This is a placeholder until metric information "
        "is available via yaml.",
        default=["count", "num", "igma", "tdev", "Repeat"],
    )
    sequentialColorMap = ListField[str](
        doc="List of hexidecimal colors for a sequential color map.",
        default=["#F5F5F5", "#5AB4AC", "#284D48"],
    )
    divergentColorMap = ListField[str](
        doc="List of hexidecimal colors for a divergent color map.",
        default=["#9A6E3A", "#C6A267", "#A9A9A9", "#4F938B", "#2C665A"],
    )

    def setDefaults(self):
        super().setDefaults()
        self.prep.keysToLoad = [
            "corners",
            "tract",
        ]
        self.produce.plot = WholeSkyPlot()

    def finalize(self):
        # Prevent finalize from running twice.
        if not self.produce.plot.plotKeys:
            # Initialize templates for assigning each key a band.
            bandSearchTemplates = [(key, key) for key in self.plotKeys]
            numBands = len(self.bands)
            # Loop over bands to define names of band-specific metrics.
            for key in self.keysWithBand:
                self.plotKeys.extend([key.format(band=band) for band in self.bands])
                # Store template and key information for band searches.
                bandSearchTemplates.extend([(key, plotKey) for plotKey in self.plotKeys[-numBands:]])
            # Assign each key a band, if any.
            keyBands = []
            for pattern, key in bandSearchTemplates:
                keyBands.append(self.findBand(pattern, key))

            # Check for duplicate keys.
            if len(set(self.plotKeys)) != len(self.plotKeys):
                raise ValueError(f"{self.plotKeys=} has duplicate keys.")
            # Check that there are actually keys to plot.
            if len(set(self.plotKeys)) == 0:
                raise ValueError(f"{self.plotKeys=} is empty.")

            # Load all config options for the plotting code.
            self.prep.keysToLoad += [key for key in self.plotKeys if key not in self.prep.keysToLoad]
            self.produce.plot.keyBands = keyBands
            self.produce.plot.plotKeys = self.plotKeys
            self.produce.plot.xAxisLabel = self.xAxisLabel
            self.produce.plot.yAxisLabel = self.yAxisLabel
            self.produce.plot.xLimits = self.xLimits
            self.produce.plot.yLimits = self.yLimits
            self.produce.plot.autoAxesLimits = self.autoAxesLimits
            self.produce.plot.figureSize = self.figureSize
            self.produce.plot.colorBarRange = self.colorBarRange
            self.produce.plot.sequentialMetrics = self.sequentialMetrics
            self.produce.plot.sequentialColorMap = self.sequentialColorMap
            self.produce.plot.divergentColorMap = self.divergentColorMap

            super().finalize()

    def validate(self):
        super().validate()
        # Check that corners and tract are in the keys to load.
        if set(("corners", "tract")) - set(self.prep.keysToLoad):
            raise ValueError(f"'corners' and 'tract' must be in {self.prep.keysToLoad=}")

    def findBand(self, pattern: str, inputString: str) -> str:
        """Search for a band in the metric key.

        Parameters
        ----------
        pattern : `str`
            A metric key with a band formatting pattern.
        inputString : `str`
            A metric key.

        Returns
        -------
        found[0] : `str`
            The band found in the search, if any.
        """
        results = [a for a in Formatter().parse(pattern)]
        finder = ""
        for i, tup in enumerate(results):
            if tup and tup[1] == "band":
                finder += tup[0]
                finder += "(.*)"
                if len(results) > i:
                    finder += results[i + 1][0]
            if not finder:
                return ""
            found = re.findall(finder, inputString)
            if found:
                return found[0]
            else:
                return ""
