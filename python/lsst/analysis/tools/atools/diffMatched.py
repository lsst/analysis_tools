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

__all__ = (
    "MatchedRefCoaddToolBase",
    "MatchedRefCoaddDiffMetric",
    "MatchedRefCoaddDiffMagTool",
    "MatchedRefCoaddDiffMagMetric",
    "MatchedRefCoaddDiffPositionTool",
    "MatchedRefCoaddDiffPositionMetric",
)

from lsst.pex.config import ChoiceField, Field

from ..actions.vector import (
    CalcBinnedStatsAction,
    ConstantValue,
    DivideVector,
    DownselectVector,
    LoadVector,
    MultiplyVector,
    SubtractVector,
)
from ..actions.vector.selectors import RangeSelector, VectorSelector
from .genericBuild import ExtendednessTool, MagnitudeXTool
from .genericProduce import MagnitudeScatterPlot


class MatchedRefCoaddToolBase(MagnitudeXTool, ExtendednessTool):
    """Base tool for matched-to-reference metrics/plots on coadds.

    Metrics/plots are expected to use the reference magnitude and
    require separate star/galaxy/all source selections.

    Notes
    -----
    The tool does not use a standard coadd flag selector, because
    it is expected that the matcher has been configured to select
    appropriate candidates (and stores a match_candidate column).
    """

    def setDefaults(self):
        super().setDefaults()
        self.mag_x = "ref_matched"

    def finalize(self):
        super().finalize()
        self._set_flux_default("mag_x")


class MatchedRefCoaddDiffTool(MatchedRefCoaddToolBase):
    """Base tool for generic diffs between reference and measured values."""

    compute_chi = Field[bool](
        default=False,
        doc="Whether to compute scaled flux residuals (chi) instead of magnitude differences",
    )

    def setDefaults(self):
        super().setDefaults()
        self.process.filterActions.yAll = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="starSelector")
        )


class MatchedRefCoaddDiffMetric(MatchedRefCoaddDiffTool):
    """Base tool for matched-to-reference metrics on coadds."""

    # These are optional because validate can be called before finalize
    # Validate should not fail in that case if it would otherwise succeed
    name_prefix = Field[str](doc="Prefix for metric key", default=None, optional=True)
    unit = Field[str](doc="Astropy unit of y-axis values", default=None, optional=True)

    _mag_low_min: int = 15
    _mag_low_max: int = 27
    _mag_interval: int = 1

    _names = ("stars", "galaxies", "all")
    _types = ("unresolved", "resolved", "all")

    def configureMetrics(
        self,
        unit: str | None = None,
        name_prefix: str | None = None,
        name_suffix: str = "_ref_mag{minimum}",
        unit_select: str = "mag",
    ):
        """Configure metric actions and return units.

        Parameters
        ----------
        unit : `str`
            The (astropy) unit of the summary statistic metrics.
        name_prefix : `str`
            The prefix for the action (column) name.
        name_suffix : `str`
            The sufffix for the action (column) name.
        unit_select : `str`
            The (astropy) unit of the selection (x-axis) column. Default "mag".

        Returns
        -------
        units : `dict` [`str`, `str`]
            A dict of the unit (value) for each metric name (key)
        """
        if unit is None:
            unit = self.unit if self.unit is not None else ""
        if name_prefix is None:
            name_prefix = self.name_prefix if self.name_prefix is not None else ""

        if unit_select is None:
            unit_select = "mag"

        key_flux = self.config_mag_x.key_flux

        units = {}
        for name, name_class in zip(self._names, self._types):
            name_capital = name.capitalize()
            x_key = f"x{name_capital}"

            for minimum in range(self._mag_low_min, self._mag_low_max + 1):
                action = getattr(self.process.calculateActions, f"{name}{minimum}")
                action.selector_range = RangeSelector(
                    vectorKey=x_key,
                    minimum=minimum,
                    maximum=minimum + self._mag_interval,
                )

                action.name_prefix = name_prefix.format(key_flux=key_flux, name_class=name_class)
                if self.parameterizedBand:
                    action.name_prefix = f"{{band}}_{action.name_prefix}"
                action.name_suffix = name_suffix.format(minimum=minimum)

                units.update(
                    {
                        action.name_median: unit,
                        action.name_sigmaMad: unit,
                        action.name_count: "count",
                        action.name_select_minimum: unit_select,
                        action.name_select_median: unit_select,
                        action.name_select_maximum: unit_select,
                    }
                )
        return units

    def setDefaults(self):
        super().setDefaults()

        for name in self._names:
            key = f"y{name.capitalize()}"
            for minimum in range(self._mag_low_min, self._mag_low_max + 1):
                setattr(
                    self.process.calculateActions,
                    f"{name}{minimum}",
                    CalcBinnedStatsAction(
                        key_vector=key,
                        selector_range=RangeSelector(
                            vectorKey=key,
                            minimum=minimum,
                            maximum=minimum + self._mag_interval,
                        ),
                    ),
                )


class MatchedRefCoaddDiffMagTool(MatchedRefCoaddDiffTool):
    """Tool for diffs between reference and measured coadd mags."""

    mag_y = Field[str](default="cmodel_err", doc="Flux (magnitude) field to difference against ref")

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            # TODO: Is this hack to ensure mag_y is set before plot tools
            # are called necessary?
            self._set_flux_default("mag_y")
            super().finalize()
            if self.compute_chi:
                self.process.buildActions.diff = DivideVector(
                    actionA=SubtractVector(
                        actionA=getattr(self.process.buildActions, f"flux_{self.mag_y}"),
                        actionB=self.process.buildActions.flux_ref_matched,
                    ),
                    actionB=getattr(self.process.buildActions, f"flux_err_{self.mag_y}"),
                )
            else:
                self.process.buildActions.diff = DivideVector(
                    actionA=SubtractVector(
                        actionA=getattr(self.process.buildActions, f"mag_{self.mag_y}"),
                        actionB=self.process.buildActions.mag_ref_matched,
                    ),
                    actionB=ConstantValue(value=1e-3),
                )


# The diamond inheritance on MatchedRefCoaddTool seems ok
class MatchedRefCoaddDiffMagMetric(MatchedRefCoaddDiffMagTool, MatchedRefCoaddDiffMetric):
    """Metric for diffs between reference and measured coadd mags."""

    def finalize(self):
        super().finalize()
        if self.unit is None:
            self.unit = "" if self.compute_chi else "mmag"
        if self.name_prefix is None:
            subtype = "chi" if self.compute_chi else "diff"
            self.name_prefix = f"photom_mag_{{key_flux}}_{{name_class}}_{subtype}_"
        if not self.produce.metric.units:
            self.produce.metric.units = self.configureMetrics()


class MatchedRefCoaddDiffMagPlot(MatchedRefCoaddDiffMagTool, MagnitudeScatterPlot):
    def finalize(self):
        # TODO: Check if this is really necessary
        # finalizing in this order should get all fluxes finalized before
        # the MagnitudeScatterPlot looks for a flux to compute S/N from
        MatchedRefCoaddDiffMagTool.finalize(self)
        MagnitudeScatterPlot.finalize(self)
        if not self.produce.yAxisLabel:
            config = self.fluxes[self.mag_y]
            label = f"{config.name_flux} - {self.fluxes['ref_matched'].name_flux}"
            self.produce.yAxisLabel = f"chi = ({label})/error" if self.compute_chi else f"{label} (mmag)"


class MatchedRefCoaddDiffPositionTool(MatchedRefCoaddDiffTool):
    """Tool for diffs between reference and measured coadd astrometry."""

    scale_factor = Field[float](
        doc="The factor to multiply positions by (i.e. the pixel scale if coordinates have pixel units)",
        default=200,
    )
    variable = ChoiceField[str](
        doc="The astrometric variable to compute metrics for",
        allowed={
            "x": "x",
            "y": "y",
        },
        optional=False,
    )

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            super().finalize()
            self.process.buildActions.pos_meas = LoadVector(vectorKey=self.variable)
            self.process.buildActions.pos_ref = LoadVector(vectorKey=f"refcat_{self.variable}")
            if self.compute_chi:
                self.process.buildActions.diff = DivideVector(
                    actionA=SubtractVector(
                        actionA=self.process.buildActions.pos_meas,
                        actionB=self.process.buildActions.pos_ref,
                    ),
                    actionB=LoadVector(vectorKey=f"{self.process.buildActions.pos_meas.vectorKey}Err"),
                )
            else:
                self.process.buildActions.diff = MultiplyVector(
                    actionA=ConstantValue(value=self.scale_factor),
                    actionB=SubtractVector(
                        actionA=self.process.buildActions.pos_meas,
                        actionB=self.process.buildActions.pos_ref,
                    ),
                )


class MatchedRefCoaddDiffPositionMetric(MatchedRefCoaddDiffPositionTool, MatchedRefCoaddDiffMetric):
    """Metric for diffs between reference and base coadd centroids."""

    def finalize(self):
        super().finalize()
        if self.unit is None:
            self.unit = "" if self.compute_chi else "mas"
        if self.name_prefix is None:
            subtype = "chi" if self.compute_chi else "diff"
            self.name_prefix = f"astrom_{self.variable}_{{name_class}}_{subtype}_"
        if not self.produce.metric.units:
            self.produce.metric.units = self.configureMetrics()


class MatchedRefCoaddDiffPositionPlot(MatchedRefCoaddDiffPositionTool, MagnitudeScatterPlot):
    # The matched catalog columns are configurable but default to cmodel only
    mag_sn = Field[str](default="cmodel_err", doc="Flux (magnitude) field to use for S/N binning on plot")

    def finalize(self):
        if not self.produce.yAxisLabel:
            # Set before MagnitudeScatterPlot.finalize or it'll default to PSF
            # Matched ref tables may not have PSF fluxes, or prefer CModel
            self._set_flux_default("mag_sn")
            super().finalize()
            self.produce.yAxisLabel = (
                f"chi = (slot - Reference {self.variable} position)/error"
                if self.compute_chi
                else f"{self.variable} position (pix)"
            )
