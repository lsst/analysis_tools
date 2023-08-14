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

from lsst.pex.config import Field, ListField

from ..actions.plot.wholeSkyPlot import WholeSkyPlot
from ..interfaces import AnalysisTool


class WholeSkyPlotTool(AnalysisTool):
    """Plot metrics across all tracts on the sky."""

    bands = ListField[str](
        doc="Photometric bands to use for plots.",
        default=["u", "g", "r", "i", "z", "y"],
        listCheck=lambda x: len(set(x)) == len(x),
    )
    parameterizedBand = Field[bool](
        doc="Does this AnalysisTool support band as a name parameter", default=False
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
    xLimits = ListField[float](doc="Plotting limits for the x axis.", default=[-5.0, 365.0])
    yLimits = ListField[float](doc="Plotting limits for the y axis.", default=[-10.0, 60.0])
    figureSize = ListField[float](doc="Size of the figure.", default=[9.0, 3.5])
    colorBarRange = Field[float](
        doc="The multiplier for the color bar range. The max/min range values are: median +/- N * sigmaMad"
        ", where N is this config value.",
        default=3.0,
    )

    propagateData = True

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
            # Loop over bands to define names of band-specific metrics.
            for key in self.keysWithBand:
                self.plotKeys.extend([key.format(band=band) for band in self.bands])
            if len(set(self.plotKeys)) != len(self.plotKeys):
                raise ValueError(f"{self.plotKeys=} has duplicate keys.")
            if len(set(self.plotKeys)) == 0:
                raise ValueError(f"{self.plotKeys=} is empty.")
            self.prep.keysToLoad += [key for key in self.plotKeys if key not in self.prep.keysToLoad]
            self.produce.plot.plotKeys = self.plotKeys
            self.produce.plot.xAxisLabel = self.xAxisLabel
            self.produce.plot.yAxisLabel = self.yAxisLabel
            self.produce.plot.xLimits = self.xLimits
            self.produce.plot.yLimits = self.yLimits
            self.produce.plot.figureSize = self.figureSize
            self.produce.plot.colorBarRange = self.colorBarRange
            super().finalize()

    def validate(self):
        if set(("corners", "tract")) - set(self.prep.keysToLoad):
            raise ValueError(f"'corners' and 'tract' must be in {self.prep.keysToLoad=}")
