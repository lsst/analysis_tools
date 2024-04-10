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

__all__ = ("MagnitudeScatterPlot",)

from lsst.pex.config import ListField

from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.vector.vectorActions import DownselectVector, VectorSelector
from .genericBuild import MagnitudeXTool


class MagnitudeScatterPlot(MagnitudeXTool):
    """A scatter plot with a magnitude on the x-axis."""

    suffixes_y_finalize = ListField[str](
        doc="Suffixes for y-axis keys that require finalization of summary stats",
        default=[""],
    )

    def setDefaults(self):
        super().setDefaults()

        # init with placeholders
        self.produce.plot = ScatterPlotWithTwoHists(xAxisLabel="", yAxisLabel="", magLabel="")
        self.produce.plot.plotTypes = ["galaxies", "stars"]
        self.produce.plot.addSummaryPlot = False

    def finalize(self):
        super().finalize()
        config_x = self.config_mag_x
        label_x = f"{config_x.name_flux} (mag)"
        # Hacky way to check if setup is complete
        if self.produce.plot.xAxisLabel == label_x:
            return
        self.produce.plot.xAxisLabel = label_x
        self.produce.plot.magLabel = self.produce.plot.xAxisLabel

        # Can't compute S/N of magnitude with no errors (e.g. true mag)
        # Try to find another or give up
        key_err = config_x.key_flux_error
        if key_err is None:
            for key, config in self.fluxes.items():
                if config.key_flux_error is not None:
                    key_err = key
                    break
            # Try to add PSF flux if all else fails
            if key_err is None:
                config_err = self.fluxes_default.psf_err
                key_err = config_err.key_flux_error
                self.fluxes["flux_sn"] = config_err
        else:
            key_err = self.mag_x

        keys_filter = [("", "flux_", self.mag_x), ("Err", "flux_err_", key_err)]
        if key_err != self.mag_x:
            keys_filter.append(("", "flux_", key_err))

        for prefix, plural in (("star", "Stars"), ("galaxy", "Galaxies")):
            for suffix, prefix_vec, key in keys_filter:
                setattr(
                    self.process.filterActions,
                    f"{prefix}_{key}_flux{suffix}",
                    DownselectVector(
                        vectorKey=f"{prefix_vec}{key}",
                        selector=VectorSelector(vectorKey=f"{prefix}Selector"),
                    ),
                )

            for suffix_y in self.suffixes_y_finalize:
                statAction = ScatterPlotStatsAction(
                    vectorKey=f"y{plural.capitalize()}{suffix_y}",
                    prefix=plural,
                    suffix=suffix_y,
                )
                fluxType = f"{prefix}_{key_err}_flux"
                statAction.highSNSelector.fluxType = fluxType
                statAction.highSNSelector.threshold = 200
                statAction.lowSNSelector.fluxType = fluxType
                statAction.lowSNSelector.threshold = 10
                statAction.fluxType = fluxType
                setattr(self.process.calculateActions, f"stats_{plural}{suffix_y}", statAction)
