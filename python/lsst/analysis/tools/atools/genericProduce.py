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
from ..actions.vector import DownselectVector
from ..actions.vector.selectors import VectorSelector
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
        # Set plot types based on config
        object_classes = []
        if self.use_galaxies:
            object_classes.append("galaxy")
        if self.use_stars:
            object_classes.append("star")
        self.produce.plot.plotTypes = [
            self.get_class_name_plural(object_class) for object_class in object_classes
        ]

        config_x = self.config_mag_x
        label_x = f"{{band}} {config_x.name_flux} (mag)"
        # Hacky way to check if setup is complete
        if self.produce.plot.xAxisLabel == label_x:
            return
        self.produce.plot.xAxisLabel = label_x
        self.produce.plot.magLabel = self.produce.plot.xAxisLabel

        # Can't compute S/N of magnitude with no errors (e.g. true mag)
        # Try to find another or give up
        key_err = config_x.key_flux_error
        name_err = config_x.name_flux_short
        if key_err is None:
            for key, config in self.fluxes.items():
                if config.key_flux_error is not None:
                    key_err = key
                    name_err = config.name_flux_short
                    break
            # Try to add PSF flux if all else fails
            if key_err is None:
                config_err = self.fluxes_default.psf_err
                key_err = config_err.key_flux_error
                name_err = config_err.name_flux_short
                self.fluxes["flux_sn"] = config_err
        else:
            key_err = self.mag_x
            name_err = self.config_mag_x.name_flux_short

        keys_filter = [
            ("", "flux_", self.mag_x, self.config_mag_x.name_flux_short),
            ("Err", "flux_err_", key_err, name_err),
        ]
        # The magnitude used for S/N is not the x-axis magnitude
        # So it has to be loaded and filtered separately
        if key_err != self.mag_x:
            keys_filter.append(("", "flux_", key_err, name_err))

        for object_class in object_classes:
            plural = self.get_class_name_plural(object_class)
            for suffix, prefix_vec, key, name_attr_mag in keys_filter:
                name_selector = self.get_name_attr_selector(object_class)
                setattr(
                    self.process.filterActions,
                    f"{object_class}_{name_attr_mag}_flux{suffix}",
                    DownselectVector(
                        vectorKey=f"{prefix_vec}{key}",
                        selector=VectorSelector(vectorKey=name_selector),
                    ),
                )

            name_y = self.get_name_attr_values(object_class)
            for suffix_y in self.suffixes_y_finalize:
                statAction = ScatterPlotStatsAction(
                    vectorKey=f"{name_y}{suffix_y}",
                    prefix=plural,
                    suffix=suffix_y,
                )
                fluxType = f"{object_class}_{name_err}_flux"
                statAction.highSNSelector.fluxType = fluxType
                statAction.highSNSelector.threshold = 200
                statAction.lowSNSelector.fluxType = fluxType
                statAction.lowSNSelector.threshold = 10
                statAction.fluxType = fluxType
                setattr(self.process.calculateActions, f"stats_{plural}{suffix_y}", statAction)
