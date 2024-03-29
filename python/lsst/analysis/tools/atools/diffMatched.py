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
    "ReferenceGalaxySelector",
    "ReferenceObjectSelector",
    "ReferenceStarSelector",
    "MatchedRefCoaddToolBase",
    "MatchedRefCoaddDiffTool",
    "MatchedRefCoaddDiffMagTool",
    "MatchedRefCoaddDiffPositionTool",
)

from lsst.pex.config import Field

from ..actions.vector import (
    CalcBinnedStatsAction,
    ConstantValue,
    DivideVector,
    DownselectVector,
    LoadVector,
    MultiplyVector,
    SubtractVector,
)
from ..actions.vector.selectors import RangeSelector, ThresholdSelector, VectorSelector
from ..interfaces import KeyedData, Vector
from .genericBuild import MagnitudeXTool
from .genericProduce import MagnitudeScatterPlot


class ReferenceGalaxySelector(ThresholdSelector):
    """A selector that selects galaxies from a catalog with a
    boolean column identifying unresolved sources.
    """

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result = super().__call__(data=data, **kwargs)
        if self.plotLabelKey:
            self._addValueToPlotInfo("true galaxies", **kwargs)
        return result

    def setDefaults(self):
        super().setDefaults()
        self.op = "eq"
        self.threshold = 0
        self.plotLabelKey = "Selection: Galaxies"
        self.vectorKey = "refcat_is_pointsource"


class ReferenceObjectSelector(RangeSelector):
    """A selector that selects all objects from a catalog with a
    boolean column identifying unresolved sources.
    """

    def setDefaults(self):
        super().setDefaults()
        self.minimum = 0
        self.vectorKey = "refcat_is_pointsource"


class ReferenceStarSelector(ThresholdSelector):
    """A selector that selects stars from a catalog with a
    boolean column identifying unresolved sources.
    """

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result = super().__call__(data=data, **kwargs)
        if self.plotLabelKey:
            self._addValueToPlotInfo("true stars", **kwargs)
        return result

    def setDefaults(self):
        super().setDefaults()
        self.op = "eq"
        self.plotLabelKey = "Selection: Stars"
        self.threshold = 1
        self.vectorKey = "refcat_is_pointsource"


class MatchedRefCoaddToolBase(MagnitudeXTool):
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
        self.process.buildActions.allSelector = ReferenceObjectSelector()
        self.process.buildActions.galaxySelector = ReferenceGalaxySelector()
        self.process.buildActions.starSelector = ReferenceStarSelector()

    def finalize(self):
        super().finalize()
        self._set_flux_default("mag_x")


class MatchedRefCoaddDiffTool(MatchedRefCoaddToolBase):
    """Base tool for generic diffs between reference and measured values."""

    compute_chi = Field[bool](
        default=False,
        doc="Whether to compute scaled flux residuals (chi) instead of magnitude differences",
    )
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

        self.process.filterActions.yAll = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="starSelector")
        )

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


class MatchedRefCoaddDiffMagTool(MatchedRefCoaddDiffTool, MagnitudeScatterPlot):
    """Tool for diffs between reference and measured coadd mags."""

    mag_y = Field[str](default="cmodel_err", doc="Flux (magnitude) field to difference against ref")

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            # Ensure mag_y is set before any plot finalizes
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
            if not self.produce.plot.yAxisLabel:
                config = self.fluxes[self.mag_y]
                label = f"{config.name_flux} - {self.fluxes['ref_matched'].name_flux}"
                self.produce.plot.yAxisLabel = (
                    f"chi = ({label})/error" if self.compute_chi else f"{label} (mmag)"
                )
            if self.unit is None:
                self.unit = "" if self.compute_chi else "mmag"
            if self.name_prefix is None:
                subtype = "chi" if self.compute_chi else "diff"
                self.name_prefix = f"photom_mag_{{key_flux}}_{{name_class}}_{subtype}_"
            if not self.produce.metric.units:
                self.produce.metric.units = self.configureMetrics()


class MatchedRefCoaddDiffPositionTool(MatchedRefCoaddDiffTool, MagnitudeScatterPlot):
    """Tool for diffs between reference and measured coadd astrometry."""

    mag_sn = Field[str](default="cmodel_err", doc="Flux (magnitude) field to use for S/N binning on plot")
    scale_factor = Field[float](
        doc="The factor to multiply positions by (i.e. the pixel scale if coordinates have pixel units)",
        default=200,
    )
    coord_label = Field[str](
        doc="The plot label for the astrometric variable (default coord_meas)",
        optional=True,
        default=None,
    )
    coord_meas = Field[str](
        doc="The key for measured values of the astrometric variable",
        optional=False,
    )
    coord_ref = Field[str](
        doc="The key for reference values of the astrometric variabler",
        optional=False,
    )

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            # Set before MagnitudeScatterPlot.finalize to undo PSF default.
            # Matched ref tables may not have PSF fluxes, or prefer CModel.
            self._set_flux_default("mag_sn")
            super().finalize()
            name = self.coord_label if self.coord_label else self.coord_meas
            self.process.buildActions.pos_meas = LoadVector(vectorKey=self.coord_meas)
            self.process.buildActions.pos_ref = LoadVector(vectorKey=self.coord_ref)
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
            if self.unit is None:
                self.unit = "" if self.compute_chi else "mas"
            if self.name_prefix is None:
                subtype = "chi" if self.compute_chi else "diff"
                self.name_prefix = f"astrom_{self.coord_meas}_{{name_class}}_{subtype}_"
            if not self.produce.metric.units:
                self.produce.metric.units = self.configureMetrics()
            if not self.produce.plot.yAxisLabel:
                self.produce.plot.yAxisLabel = (
                    f"chi = (slot - reference {name} position)/error"
                    if self.compute_chi
                    else f"slot - reference {name} position ({self.unit})"
                )
