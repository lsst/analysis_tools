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

__all__ = ("CompletenessPerPatchPropertyMapPlot",)

from ..actions.plot.patchActionSkyPlot import PerPatchMetricConfig, PerPatchPropertyMapPlot


class CompletenessPerPatchPropertyMapPlot(PerPatchPropertyMapPlot):
    """Per-patch completeness plots from a matched catalog."""

    def setDefaults(self):
        self.cmap_metric = "RdYlBu"
        self.metrics = {
            "mag_compl50_any": PerPatchMetricConfig(
                description="Magnitude at 50% completeness",
                key="{band}_detect_{name_flux_target}_vs_{name_flux_ref}_all_mag_completeness_50p00_pct",
                vmin=25.0,
                vmax=27.5,
            ),
            "mag_compl80_any": PerPatchMetricConfig(
                description="Magnitude at 80% completeness",
                key="{band}_detect_{name_flux_target}_vs_{name_flux_ref}_all_mag_completeness_80p00_pct",
                vmin=25,
                vmax=26,
            ),
            "mag_compl90_any": PerPatchMetricConfig(
                description="Magnitude at 90% completeness",
                key="{band}_detect_{name_flux_target}_vs_{name_flux_ref}_all_mag_completeness_90p00_pct",
                vmin=24,
                vmax=25,
            ),
            "compl_24_25_any": PerPatchMetricConfig(
                description="Completeness for 24 < reference mag < 25 objects",
                key="{band}_detect_{name_flux_target}_vs_{name_flux_ref}_all_completeness_mag24p0",
                vmin=0.80,
                vmax=0.95,
            ),
        }
