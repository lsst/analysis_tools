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
    "MatchedRefCoaddDiffColorTool",
    "MatchedRefCoaddDiffMagTool",
    "MatchedRefCoaddDiffPositionTool",
)

import copy
from abc import abstractmethod

from lsst.pex.config import DictField, Field

from ..actions.vector import (
    CalcBinnedStatsAction,
    ColorDiff,
    ColorError,
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
from .genericMetricAction import StructMetricAction
from .genericPlotAction import StructPlotAction
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

    _names = {"stars": "star", "galaxies": "galaxy", "all": "all"}
    _types = ("unresolved", "resolved", "all")

    def _set_actions(self, suffix=None):
        if suffix is None:
            suffix = ""

        for name_lower, name_singular in self._names.items():
            name = name_lower.capitalize()
            key = f"y{name}{suffix}"
            setattr(
                self.process.filterActions,
                key,
                DownselectVector(
                    vectorKey=f"diff{suffix}", selector=VectorSelector(vectorKey=f"{name_singular}Selector")
                ),
            )

            for minimum in range(self._mag_low_min, self._mag_low_max + 1):
                setattr(
                    self.process.calculateActions,
                    f"{name_lower}{minimum}{suffix}",
                    CalcBinnedStatsAction(
                        key_vector=key,
                        selector_range=RangeSelector(
                            vectorKey=key,
                            minimum=minimum,
                            maximum=minimum + self._mag_interval,
                        ),
                        return_minmax=False,
                    ),
                )

    def configureMetrics(
        self,
        unit: str | None = None,
        name_prefix: str | None = None,
        name_suffix: str = "_ref_mag{minimum}",
        attr_suffix: str | None = None,
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
            The suffix for the action (column) name.
        attr_suffix : `str`
            The suffix for the attribute to assign the action to.
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
        if attr_suffix is None:
            attr_suffix = ""

        if unit_select is None:
            unit_select = "mag"

        key_flux = self.config_mag_x.key_flux

        units = {}
        for name, name_class in zip(self._names, self._types):
            name_capital = name.capitalize()
            x_key = f"x{name_capital}"

            for minimum in range(self._mag_low_min, self._mag_low_max + 1):
                action = getattr(self.process.calculateActions, f"{name}{minimum}{attr_suffix}")
                action.selector_range = RangeSelector(
                    vectorKey=x_key,
                    minimum=minimum,
                    maximum=minimum + self._mag_interval,
                )

                action.name_prefix = name_prefix.format(
                    key_flux=key_flux,
                    name_class=name_class,
                )
                if self.parameterizedBand:
                    action.name_prefix = f"{{band}}_{action.name_prefix}"
                action.name_suffix = name_suffix.format(minimum=minimum)

                units.update(
                    {
                        action.name_median: unit,
                        action.name_sigmaMad: unit,
                        action.name_count: "count",
                        action.name_select_median: unit_select,
                    }
                )
        return units

    @property
    def config_mag_y(self):
        """Return the y-axis magnitude config.

        Although tools may not end up using any flux measures in metrics or
        plots, this should still be set to the flux measure that was matched
        or selected against in the catalog not used for the x-axis."""
        mag_y = self.get_key_flux_y()
        if mag_y not in self.fluxes:
            raise KeyError(f"{mag_y=} not in {self.fluxes}; was finalize called?")
        # This is a logic error: it shouldn't be called before finalize
        assert mag_y in self.fluxes
        return self.fluxes[mag_y]

    @abstractmethod
    def get_key_flux_y(self) -> str:
        """Return the key for the y-axis flux measure."""
        raise NotImplementedError("subclasses must implement get_key_flux_y")

    def setDefaults(self):
        super().setDefaults()
        self._set_actions()


class MatchedRefCoaddDiffColorTool(MatchedRefCoaddDiffTool, MagnitudeScatterPlot):
    """Tool for diffs between reference and measured coadd mags.

    Notes
    -----
    Since this tool requires at least two bands, it is essentially impossible
    to call on its own.
    """

    mag_y1 = Field[str](default="cmodel_err", doc="Flux field for first magnitude")
    mag_y2 = Field[str](
        doc="Flux field for second magnitude (to subtract from first); default same as first",
        default=None,
        optional=True,
    )
    bands = DictField[str, str](
        doc="Bands for colors. ",
        default={"u": "g", "g": "r,i", "r": "i", "i": "z", "z": "y"},
    )
    band_separator = Field[str](default=",", doc="Separator for multiple bands")

    def _split_bands(self, band_list: str):
        return band_list.split(self.band_separator)

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff_0"):
            if self.mag_y2 is None:
                self.mag_y2 = self.mag_y1
            # Ensure mag_y1/2 are set before any plot finalizes
            # This may result in duplicate actions but these are just plain
            # column selectors so that's not a serious problem
            bands = {band1: self._split_bands(band2_list) for band1, band2_list in self.bands.items()}
            n_bands = 0

            for band1, band2_list in bands.items():
                for band2 in band2_list:
                    mag_y1 = f"mag_y_{band1}"
                    mag_y2 = f"mag_y_{band2}"
                    mag_x1 = f"mag_x_{band1}"
                    mag_x2 = f"mag_x_{band2}"
                    self._set_flux_default(mag_y1, band=band1, name_mag=self.mag_y1)
                    self._set_flux_default(mag_y2, band=band2, name_mag=self.mag_y2)
                    self._set_flux_default(mag_x1, band=band1, name_mag=self.mag_x)
                    self._set_flux_default(mag_x2, band=band2, name_mag=self.mag_x)
                    n_bands += 1

            self.suffixes_y_finalize = [f"_{idx}" for idx in range(n_bands)]
            super().finalize()

            self.unit = "" if self.compute_chi else "mmag"
            subtype = "chi" if self.compute_chi else "diff"

            metric_base = self.produce.metric
            plot_base = self.produce.plot

            actions_metric = {}
            actions_plot = {}

            config_mag_x = self.config_mag_x
            config_mag_y = self.config_mag_y
            name_short_x = config_mag_x.name_flux_short
            name_short_y = config_mag_y.name_flux_short

            idx = 0
            for band1, band2_list in bands.items():
                for band2 in band2_list:
                    name_color = f"{band1}_minus_{band2}"
                    # Keep this index-based to simplify finalize
                    suffix_y = f"_{idx}"
                    self._set_actions(suffix=suffix_y)
                    metric = copy.copy(metric_base)
                    self.name_prefix = (
                        f"photom_{name_short_y}_vs_{name_short_x}_color_{name_color}"
                        f"_{subtype}_{{name_class}}_"
                    )
                    metric.units = self.configureMetrics(attr_suffix=suffix_y)
                    plot = copy.copy(plot_base)

                    plot.suffix_y = suffix_y
                    plot.suffix_stat = suffix_y

                    mag_y1 = f"{self.mag_y1}_{band1}"
                    mag_y2 = f"{self.mag_y2}_{band2}"
                    mag_x1 = f"{self.mag_x}_{band1}"
                    mag_x2 = f"{self.mag_x}_{band2}"

                    diff = ColorDiff(
                        color1_flux1=getattr(self.process.buildActions, f"flux_{mag_y1}"),
                        color1_flux2=getattr(self.process.buildActions, f"flux_{mag_y2}"),
                        color2_flux1=getattr(self.process.buildActions, f"flux_{mag_x1}"),
                        color2_flux2=getattr(self.process.buildActions, f"flux_{mag_x2}"),
                    )

                    if self.compute_chi:
                        diff = DivideVector(
                            actionA=diff,
                            actionB=ColorError(
                                flux_err1=DivideVector(
                                    actionA=getattr(self.process.buildActions, f"flux_err_{mag_y1}"),
                                    actionB=getattr(self.process.buildActions, f"flux_{mag_y1}"),
                                ),
                                flux_err2=DivideVector(
                                    actionA=getattr(self.process.buildActions, f"flux_err_{mag_y2}"),
                                    actionB=getattr(self.process.buildActions, f"flux_{mag_y2}"),
                                ),
                            ),
                        )
                    setattr(self.process.buildActions, f"diff{plot.suffix_y}", diff)

                    label = f"({band1} - {band2}) ({config_mag_y.name_flux} - {config_mag_x.name_flux})"
                    label = f"chi = ({label})/error" if self.compute_chi else f"{label} (mmag)"
                    plot.yAxisLabel = label
                    actions_metric[name_color] = metric
                    actions_plot[name_color] = plot
                    idx += 1
            action_metric = StructMetricAction()
            for name_action, action in actions_metric.items():
                setattr(action_metric.actions, name_action, action)
            self.produce.metric = action_metric
            action_plot = StructPlotAction()
            for name_action, action in actions_plot.items():
                setattr(action_plot.actions, name_action, action)
            self.produce.plot = action_plot

    def get_key_flux_y(self) -> str:
        return self.mag_y1

    def setDefaults(self):
        # skip MatchedRefCoaddDiffTool.setDefaults's _setActions call
        MatchedRefCoaddToolBase.setDefaults(self)
        MagnitudeScatterPlot.setDefaults(self)

    def validate(self):
        super().validate()
        errors = []
        for band1, band2_list in self.bands.items():
            bands = self._split_bands(band2_list)
            if len(set(bands)) != len(bands):
                errors.append(f"value={band2_list} is not a set for key={band1}")
        if errors:
            raise ValueError("\n".join(errors))


class MatchedRefCoaddDiffMagTool(MatchedRefCoaddDiffTool, MagnitudeScatterPlot):
    """Tool for diffs between reference and measured coadd mags."""

    mag_y = Field[str](default="cmodel_err", doc="Flux (magnitude) field to difference against ref")

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            # Ensure mag_y is set before any plot finalizes
            self._set_flux_default("mag_y")
            super().finalize()
            name_short_x = self.config_mag_x.name_flux_short
            name_short_y = self.config_mag_y.name_flux_short

            prefix_action = "flux" if self.compute_chi else "mag"
            action_diff = SubtractVector(
                actionA=getattr(self.process.buildActions, f"{prefix_action}_{self.mag_x}"),
                actionB=getattr(self.process.buildActions, f"{prefix_action}_{self.mag_y}"),
            )

            if self.compute_chi:
                key_err = f"flux_err_{self.mag_y}"
                action_err = (
                    getattr(self.process.buildActions, key_err)
                    if hasattr(self.process.buildActions, key_err)
                    else getattr(self.process.buildActions, f"flux_err_{self.mag_x}")
                )
                self.process.buildActions.diff = DivideVector(actionA=action_diff, actionB=action_err)
            else:
                # set to mmag
                self.process.buildActions.diff = MultiplyVector(
                    actionA=action_diff,
                    actionB=ConstantValue(value=1000.0),
                )
            if not self.produce.plot.yAxisLabel:
                label = f"{self.config_mag_y.name_flux} - {self.config_mag_x.name_flux}"
                self.produce.plot.yAxisLabel = (
                    f"chi = ({label})/error" if self.compute_chi else f"{label} (mmag)"
                )
            if self.unit is None:
                self.unit = "" if self.compute_chi else "mmag"
            if self.name_prefix is None:
                subtype = "chi" if self.compute_chi else "diff"
                self.name_prefix = f"photom_{name_short_y}_vs_{name_short_x}_mag_{subtype}_{{name_class}}_"
            if not self.produce.metric.units:
                self.produce.metric.units = self.configureMetrics()

    def get_key_flux_y(self) -> str:
        return self.mag_y


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
            name_short_x = self.config_mag_x.name_flux_short
            name_short_y = self.config_mag_y.name_flux_short

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
                self.name_prefix = (
                    f"astrom_{name_short_y}_vs_{name_short_x}_{self.coord_meas}_coord_{subtype}"
                    f"_{{name_class}}_"
                )
            if not self.produce.metric.units:
                self.produce.metric.units = self.configureMetrics()
            if not self.produce.plot.yAxisLabel:
                label = f"({name_short_y} - {name_short_x})"
                self.produce.plot.yAxisLabel = (
                    f"chi = ({label} {name} coord)/error"
                    if self.compute_chi
                    else f"{label} {name} coord ({self.unit})"
                )

    def get_key_flux_y(self) -> str:
        return self.mag_sn
