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
    "MatchedRefCoaddTool",
    "MatchedRefCoaddChiColorTool",
    "MatchedRefCoaddChiCoordDecTool",
    "MatchedRefCoaddChiCoordRaTool",
    "MatchedRefCoaddChiDistanceTool",
    "MatchedRefCoaddChiMagTool",
    "MatchedRefCoaddCompurityTool",
    "MatchedRefCoaddDiffColorTool",
    "MatchedRefCoaddDiffColorZoomTool",
    "MatchedRefCoaddDiffCoordDecTool",
    "MatchedRefCoaddDiffCoordDecZoomTool",
    "MatchedRefCoaddDiffCoordRaTool",
    "MatchedRefCoaddDiffCoordRaZoomTool",
    "MatchedRefCoaddDiffDistanceTool",
    "MatchedRefCoaddDiffMagTool",
    "MatchedRefCoaddDiffMagZoomTool",
    "MatchedRefCoaddDiffPositionTool",
    "MatchedRefCoaddDiffTool",
    "MatchedRefCoaddDiffDistanceZoomTool",
    "reconfigure_diff_matched_defaults",
)

import copy
import inspect
from abc import abstractmethod

import astropy.units as u
import lsst.pex.config as pexConfig
from lsst.pex.config import DictField, Field
from lsst.pex.config.configurableActions import ConfigurableActionField

from ..actions.config import MagnitudeBinConfig
from ..actions.keyedData import (
    CalcBinnedCompletenessAction,
    CalcCompletenessHistogramAction,
    MagnitudeCompletenessConfig,
)
from ..actions.plot import CompletenessHist
from ..actions.vector import (
    CalcBinnedStatsAction,
    ColorDiff,
    ColorError,
    ConstantValue,
    CosVector,
    DivideVector,
    DownselectVector,
    IsMatchedObjectSameClass,
    LoadVector,
    MultiplyVector,
    SubtractVector,
)
from ..actions.vector.selectors import (
    InjectedGalaxySelector,
    InjectedObjectSelector,
    InjectedStarSelector,
    MatchedObjectSelector,
    RangeSelector,
    ReferenceGalaxySelector,
    ReferenceObjectSelector,
    ReferenceStarSelector,
    SelectorBase,
    VectorSelector,
)
from ..interfaces import AnalysisBaseConfig, BaseMetricAction, NoMetric
from .genericBuild import MagnitudeTool, MagnitudeXTool, ObjectClassTool
from .genericMetricAction import StructMetricAction
from .genericPlotAction import StructPlotAction
from .genericProduce import MagnitudeScatterPlot


class MatchedRefCoaddTool(ObjectClassTool):
    """Base tool for matched-to-reference metrics/plots on coadds.

    This tool is designed to configure plots and metrics as a function of
    magnitude (object or reference). The metrics are binned by the same
    magnitude shown on the x-axis in plots. By default, this is the reference
    magnitude but plots can be configured to bin by object magnitude instead.

    Notes
    -----
    The tool does not use a standard coadd flag selector, because
    it is expected that the matcher has been configured to select
    appropriate candidates (and stores a match_candidate column).

    The tool requires specification of reference galaxy and star selectors,
    as these will be used to determine whether matched objects have the same
    class as the reference, even if a particular class is not being plotted.
    It is okay to specify a "dummy" selector that always returns False if
    there are no reference objects of the given class.
    """

    _suffix_ref = "_ref"
    _suffix_target = "_target"

    context = pexConfig.ChoiceField[str](
        doc="The context for the selectors",
        allowed={
            "custom": "User-configured selectors",
            "DC2": "DC2 Truth Summary match",
            "injection": "Source injection match",
        },
        default="DC2",
    )

    select_ref_by_default = pexConfig.Field[bool](
        doc="Whether reference quantities should be used by default in other tools,"
        " e.g. for binning metrics and for the x-axis in plots",
        default=True,
    )

    selector_ref_all = ConfigurableActionField[SelectorBase](
        doc="The selector for reference objects of all types",
        default=ReferenceObjectSelector,
    )
    selector_ref_galaxy = ConfigurableActionField[SelectorBase](
        doc="The selector for reference galaxies",
        default=ReferenceGalaxySelector,
    )
    selector_ref_star = ConfigurableActionField[SelectorBase](
        doc="The selector for reference stars",
        default=ReferenceStarSelector,
    )

    mag_bins = pexConfig.ConfigField[MagnitudeBinConfig](doc="Magnitude bin configuration for metrics")
    # These are optional because validate can be called before finalize
    # Validate should not fail in that case if it would otherwise succeed
    name_prefix = pexConfig.Field[str](
        doc="Default prefix for metric key. Can include {name_type} as a"
        " template for the type of object (resolved/unresolved)",
        default=None,
        optional=True,
    )
    name_suffix = pexConfig.Field[str](
        doc="The suffix for metric names. Can include {name_mag} as a "
        " template for the magnitude algorithm",
        default="_ref_mag{name_mag}",
    )
    unit = pexConfig.Field[str](doc="Astropy unit of y-axis values", default=None, optional=True)

    def finalize(self):
        # Don't do anything if the value is the one for which the defaults of
        # selector_ref_all, etc are - this can't easily be inferred and must
        # be kept in sync manually
        if self.context != "DC2":
            match self.context:
                case "injection":
                    self.selector_ref_all = InjectedObjectSelector()
                    self.selector_ref_galaxy = InjectedGalaxySelector()
                    self.selector_ref_star = InjectedStarSelector()
                case "custom":
                    pass
                case _:
                    raise NotImplementedError(f"{self.context=} is not implemented in {self.__class__}")

        # Other tools will except selector_all
        self.selection_suffix = self._suffix_ref if self.select_ref_by_default else self._suffix_target

        super().finalize()

        for object_class in self.get_classes():
            name_selector = self.get_name_attr_selector(object_class, self._suffix_ref)
            selector = self.get_selector_ref(object_class)
            # This is a build action because selectors in prep are applied with
            # and; we're not using these to filter all points but to make
            # several parallel selections
            setattr(self.process.buildActions, name_selector, selector)

    def get_selector_ref(self, object_class: str):
        match object_class:
            case "any":
                return self.selector_ref_all
            case "galaxy":
                return self.selector_ref_galaxy
            case "star":
                return self.selector_ref_star

    def reconfigure(
        self,
        context: str | None = None,
        key_flux_meas: str | None = None,
        bands_color: dict[str, str] | list[str] | None = None,
        use_any: bool | None = None,
        use_galaxies: bool | None = None,
        use_stars: bool | None = None,
    ):
        """Reconfigure any MatchedRefCoaddTools in an analysis task config.

        Parameters
        ----------
        context
            The context to set. Must be a valid choice for
            MatchedRefCoaddTool.context.
        key_flux_meas
            The key of the measured flux config to use, e.g. "psf". If the key
            is not found, it will search for f"{key}_err", the default name for
            configurations that load error keys as well as fluxes.
        bands_color
            A dictionary keyed by band of comma-separated bands to measure
            colors for, where the color is (key - value). If a list is passed,
            tools will modify the defaults to select only those bands within
            the list (which should also be a set).
        use_any
            Whether to compute metrics for objects of all types.
        use_galaxies
            Whether to compute metrics and plot lines for galaxies only.
        use_stars
            Whether to compute metrics and plot lines for stars only.

        Notes
        -----
        Any kwargs set to None will not change the relevant config fields.
        """

        if context is not None:
            self.context = context
        if use_any is not None:
            self.use_any = use_any
        if use_galaxies is not None:
            self.use_galaxies = use_galaxies
        if use_stars is not None:
            self.use_stars = use_stars

        # This allows the method to work automatically on class defaults
        kwargs = {"self": self} if inspect.isclass(self) else {}

        # Change any dependent magnitudes
        self.reconfigure_dependent_magnitudes(key_flux_meas=key_flux_meas, bands_color=bands_color, **kwargs)

    def reconfigure_dependent_magnitudes(
        self,
        key_flux_meas: str | None = None,
        bands_color: dict[str, str] | list[str] | None = None,
    ):
        """Reconfigure any dependent (i.e., on the y-axis in plots) magnitude
           column configs.

        Parameters
        ----------
        key_flux_meas
            The key of the measured flux config to set to.
        bands_color
            A dictionary keyed by band of comma-separated bands to measure
            colors for, where the color is (key - value). If a list is passed,
            tools will modify the defaults to select only those bands within
            the list (which should also be a set).
        """

    def setDefaults(self):
        super().setDefaults()
        # The selection info isn't useful in plots with multiple classes
        self.selector_ref_galaxy.plotLabelKey = None
        self.selector_ref_star.plotLabelKey = None


class MatchedRefCoaddDiffTool(MagnitudeXTool, MatchedRefCoaddTool):
    """Base tool for generic diffs between reference and measured values."""

    limits_chi_default = (-5, 5)
    limits_diff_color_mmag_default = (-250.0, 250.0)
    limits_diff_color_mmag_zoom_default = (-50.0, 50.0)
    limits_diff_mag_mmag_default = (-1000.0, 1000.0)
    limits_diff_mag_mmag_zoom_default = (-50.0, 50.0)
    limits_diff_pos_mas_default = (-500, 500)
    limits_diff_pos_mas_zoom_default = (-10, 10)
    limits_x_mag_default = (16.0, 31.0)
    limits_x_mag_zoom_default = (16.0, 23.5)

    compute_chi = pexConfig.Field[bool](
        default=False,
        doc="Whether to compute scaled flux residuals (chi) instead of magnitude differences",
    )

    def _set_actions(self, suffix=None):
        if suffix is None:
            suffix = ""

        selection = self._suffix_ref if self.select_ref_by_default else self._suffix_target
        for object_class in self.get_classes():
            name_type_plural = self.get_class_name_plural(object_class)
            name_attr = f"{self.get_name_attr_values(object_class)}{suffix}"
            name_selector = self.get_name_attr_selector(object_class, selection)
            name_x = f"x{name_type_plural.capitalize()}"

            y_values = DownselectVector(
                vectorKey=f"diff{suffix}",
                selector=VectorSelector(vectorKey=name_selector),
            )
            setattr(self.process.filterActions, name_attr, y_values)

            bins = self.mag_bins.get_bins()
            for minimum in bins:
                setattr(
                    self.process.calculateActions,
                    f"{name_type_plural}_{minimum}{suffix}",
                    CalcBinnedStatsAction(
                        key_vector=name_attr,
                        selector_range=RangeSelector(
                            vectorKey=name_x,
                            minimum=minimum,
                            maximum=minimum + self.mag_bins.mag_width,
                        ),
                    ),
                )

    def configureMetrics(
        self,
        unit: str | None = None,
        name_prefix: str | None = None,
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

        for object_class in self.get_classes():
            name_type = self.get_class_type(object_class)
            name_type_plural = self.get_class_name_plural(object_class)
            name_capital = name_type_plural.capitalize()
            x_key = f"x{name_capital}"

            # Set up metrics for objects of one class within a magnitude range
            bins = self.mag_bins.get_bins()
            for minimum in bins:
                action = getattr(self.process.calculateActions, f"{name_type_plural}_{minimum}{attr_suffix}")
                action.selector_range = RangeSelector(
                    vectorKey=x_key,
                    minimum=minimum / 1000.0,
                    maximum=(minimum + self.mag_bins.mag_width) / 1000.0,
                )
                name_mag = self.mag_bins.get_name_bin(minimum)

                action.name_prefix = name_prefix.format(
                    key_flux=key_flux,
                    name_type=name_type,
                )
                if self.parameterizedBand:
                    action.name_prefix = f"{{band}}_{action.name_prefix}"
                action.name_suffix = self.name_suffix.format(name_mag=name_mag)

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

    def finalize(self):
        MagnitudeXTool.finalize(self)
        MatchedRefCoaddTool.finalize(self)

    @abstractmethod
    def get_key_flux_y(self) -> str:
        """Return the key for the y-axis flux measure."""
        raise NotImplementedError("subclasses must implement get_key_flux_y")

    def setDefaults(self):
        MagnitudeXTool.setDefaults(self)
        MatchedRefCoaddTool.setDefaults(self)
        self.mag_x = "ref_matched"
        self.prep.selectors.matched = MatchedObjectSelector()


class MatchedRefCoaddDiffPlot(MatchedRefCoaddDiffTool, MagnitudeScatterPlot):
    """Base tool for generic diffs between reference and measured values,
    with a scatter plot."""

    def do_metrics(self):
        return not isinstance(self.produce.metric, NoMetric)

    def get_key_flux_y(self) -> str:
        return super().get_key_flux_y()

    def finalize(self):
        MatchedRefCoaddDiffTool.finalize(self)
        MagnitudeScatterPlot.finalize(self)

    def setDefaults(self):
        # This will set no plot
        MatchedRefCoaddDiffTool.setDefaults(self)
        # This will set the plot
        MagnitudeScatterPlot.setDefaults(self)
        self.produce.plot.xLims = self.limits_x_mag_default


class MatchedRefCoaddCompurityTool(MagnitudeTool, MatchedRefCoaddTool):
    """Plot the fraction of injected sources recovered by input magnitude.

    By contrast with MatchedRefCoaddDiffTool, where one must choose which
    magnitude appears on the x-axis, this tools creates two plots with
    different magnitudes. The completeness plot necessarily is a function
    of reference magnitude while purity is a function of object (target)
    magnitude.
    """

    config_metrics = pexConfig.ConfigField[MagnitudeCompletenessConfig](
        doc="Plot-based (unbinned) metric definition configuration"
    )
    key_match_distance = pexConfig.Field[str](
        default="match_distance",
        doc="Key for match distance column (>=0 for a successful match)",
    )
    mag_bins_plot = pexConfig.ConfigField[MagnitudeBinConfig](
        doc="Magnitude bin configuration for plots and for unbinned metrics"
        "(including completeness at magnitude thresholds)"
    )
    mag_ref = pexConfig.Field[str](
        default="ref_matched",
        doc="Flux (magnitude) config key  (to self.fluxes) for reference (true) magnitudes",
    )
    mag_target = pexConfig.Field[str](
        default="cmodel_err",
        doc="Flux (magnitude) config key (to self.fluxes) for target (measured) magnitudes",
    )
    make_plots = pexConfig.Field[bool](
        default=True,
        doc="Whether to generate plots in addition to metrics",
    )

    @property
    def config_mag_ref(self):
        return self._config_mag("mag_ref")

    @property
    def config_mag_target(self):
        return self._config_mag("mag_target")

    def finalize(self):
        if not self.produce.metric.units:
            MagnitudeTool.finalize(self)
            MatchedRefCoaddTool.finalize(self)
            self._set_flux_default("mag_ref")
            self._set_flux_default("mag_target")

            # This is the default convention for metric names, originally set
            # for DC2 truth match but expanded to generic reference catalogs
            # (including injection catalogs)
            name_prefix = (
                self.name_prefix
                if self.name_prefix
                else (
                    f"detect_{self.config_mag_target.name_flux_short}_vs_"
                    f"{self.config_mag_ref.name_flux_short}_{{name_type}}_"
                )
            )
            unit_select = ""
            kwargs_matched_class_action = {}

            # Set up selectors for all object classes as they may be needed by
            # the wrong/right matched class selector
            for object_class in ("any", "galaxy", "star"):
                for suffix, func_selector in (
                    (self._suffix_ref, self.get_selector_ref),
                    (self._suffix_target, self.get_selector),
                ):
                    name_selector = self.get_name_attr_selector(object_class, suffix)
                    if not hasattr(self.process.buildActions, name_selector):
                        selector = func_selector(object_class)
                        setattr(self.process.buildActions, name_selector, selector)
                    if object_class != "any":
                        kwargs_matched_class_action[f"key_is{suffix}_{object_class}"] = name_selector

            # This isn't exactly a filterAction but by default it needs to go
            # after build and before calc, so here it is
            self.process.filterActions.matched_class = IsMatchedObjectSameClass(**kwargs_matched_class_action)

            key_flux = self.config_mag_ref.key_flux
            key_mag_ref = f"mag_{self.mag_ref}"
            key_mag_target = f"mag_{self.mag_target}"
            object_classes = self.get_classes()
            self.produce.metric = StructMetricAction()
            if self.make_plots:
                self.produce.plot = StructPlotAction()

            for object_class in object_classes:
                name_type = self.get_class_type(object_class)
                name_selector_ref = self.get_name_attr_selector(object_class, self._suffix_ref)
                name_selector_target = self.get_name_attr_selector(object_class, self._suffix_target)
                name_prefix_class = name_prefix.format(
                    key_flux=key_flux,
                    name_type=name_type,
                )
                if self.parameterizedBand:
                    name_prefix_class = f"{{band}}_{name_prefix_class}"

                units = {}
                completeness_binned_metrics = CalcCompletenessHistogramAction(
                    action=CalcBinnedCompletenessAction(
                        name_prefix=name_prefix_class,
                        selector_range_ref=RangeSelector(vectorKey=key_mag_ref),
                        selector_range_target=RangeSelector(vectorKey=key_mag_target),
                        key_mask_ref=name_selector_ref,
                        key_mask_target=name_selector_target,
                    ),
                    bins=self.mag_bins,
                )
                # Metric bins should be coarser than plot bins and therefore
                # are unsuited for computing unbinned metrics (like mag at a
                # given completeness/purity)
                completeness_binned_metrics.config_metrics.completeness_percentiles = []
                setattr(
                    self.process.calculateActions,
                    f"completeness_binned_metrics_{object_class}",
                    completeness_binned_metrics,
                )

                bins = self.mag_bins.get_bins()
                for minimum in bins:
                    name_mag = self.mag_bins.get_name_bin(minimum)
                    action = CalcBinnedCompletenessAction(
                        name_prefix=name_prefix_class,
                        name_suffix=self.name_suffix.format(name_mag=name_mag),
                        selector_range_ref=RangeSelector(
                            vectorKey=key_mag_ref,
                            minimum=minimum / 1000.0,
                            maximum=(minimum + self.mag_bins.mag_width) / 1000.0,
                        ),
                        selector_range_target=RangeSelector(
                            vectorKey=key_mag_target,
                            minimum=minimum / 1000.0,
                            maximum=(minimum + self.mag_bins.mag_width) / 1000.0,
                        ),
                        key_mask_ref=name_selector_ref,
                        key_mask_target=name_selector_target,
                    )
                    setattr(
                        self.process.calculateActions,
                        f"completeness_{object_class}_{minimum}",
                        action,
                    )

                    units.update(
                        {
                            action.name_count: "count",
                            action.name_count_ref: "count",
                            action.name_count_target: "count",
                            action.name_completeness: unit_select,
                            action.name_completeness_bad_match: unit_select,
                            action.name_completeness_good_match: unit_select,
                            action.name_purity: unit_select,
                            action.name_purity_bad_match: unit_select,
                            action.name_purity_good_match: unit_select,
                        }
                    )

                completeness_plot = CalcCompletenessHistogramAction(
                    action=CalcBinnedCompletenessAction(
                        name_prefix=name_prefix_class,
                        selector_range_ref=RangeSelector(vectorKey=key_mag_ref),
                        selector_range_target=RangeSelector(vectorKey=key_mag_target),
                        key_mask_ref=name_selector_ref,
                        key_mask_target=name_selector_target,
                    ),
                    bins=self.mag_bins_plot,
                    config_metrics=self.config_metrics,
                )
                setattr(
                    self.process.calculateActions,
                    f"completeness_plot_{object_class}",
                    completeness_plot,
                )
                for pct in completeness_plot.config_metrics.completeness_percentiles:
                    name_pct = completeness_plot.action.name_mag_completeness(
                        completeness_plot.getPercentileName(pct)
                    )
                    units[name_pct] = unit_select

                # Make the metric action for the given object class
                # This will include units for metrics from the plot histogram
                # (i.e. the magnitude for a given completeness threshold)
                setattr(
                    self.produce.metric.actions,
                    object_class,
                    BaseMetricAction(units=units),
                )

                if self.make_plots:
                    setattr(
                        self.produce.plot.actions,
                        object_class,
                        CompletenessHist(action=completeness_plot),
                    )

    def reconfigure_dependent_magnitudes(
        self,
        key_flux_meas: str | None = None,
        bands_color: dict[str, str] | list[str] | None = None,
    ):
        if key_flux_meas is not None:
            self.mag_target = key_flux_meas

    def setDefaults(self):
        MagnitudeTool.setDefaults(self)
        MatchedRefCoaddTool.setDefaults(self)

        self.mag_bins_plot.mag_interval = 100
        self.mag_bins_plot.mag_width = 200
        # Completeness/purity don't need a ref/target suffix as they are by
        # definition a function of ref/target mags, respectively
        self.name_suffix = "_mag{name_mag}"


class MatchedRefCoaddDiffColorTool(MatchedRefCoaddDiffPlot):
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
        # The empty value for y is needed to indicate that it's a valid band
        default={"u": "g", "g": "r,i", "r": "i", "i": "z", "z": "y", "y": ""},
    )
    band_separator = Field[str](default=",", doc="Separator for multiple bands")

    def _split_bands(self, band_list: str):
        # Split returns [""] for an empty string
        return band_list.split(self.band_separator) if band_list else []

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

            # Set up mag actions for every band needed before finalizing plots
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

            # These two lines must appear in this order so that every color
            # has its plot actions finalized with a suffix (i.e., pointing
            # summary stats at yStars_0 instead of yStars).
            self.suffixes_y_finalize = [f"_{idx}" for idx in range(n_bands)]
            super().finalize()

            self.unit = "" if self.compute_chi else "mmag"
            subtype = "chi" if self.compute_chi else "diff"

            metric_base = self.produce.metric
            metric = metric_base
            plot_base = self.produce.plot

            do_metrics = self.do_metrics()

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
                    self.name_prefix = (
                        f"photom_{name_short_y}_vs_{name_short_x}_color_{name_color}"
                        f"_{subtype}_{{name_type}}_"
                    )
                    if do_metrics:
                        metric = copy.copy(metric_base)
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
            if do_metrics:
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

    def reconfigure_dependent_magnitudes(
        self,
        key_flux_meas: str | None = None,
        bands_color: dict[str, str] | list[str] | None = None,
    ):
        if key_flux_meas is not None:
            self.mag_y1 = key_flux_meas
        if bands_color is not None:
            if isinstance(bands_color, dict):
                self.bands = bands_color
            else:
                bands_new = {}
                bands_old = self.bands.default if inspect.isclass(self) else self.bands
                for band in bands_color:
                    colors = bands_old.get(band)
                    if colors is None:
                        raise ValueError(
                            f"Passed {bands_color=} to reconfigure colors for {self=} but {band=}"
                            f" is not in {bands_old=}."
                        )
                    bands_new[band] = ",".join(band for band in colors.split(",") if band in bands_color)
                self.bands = bands_new

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_color_mmag_default

    def validate(self):
        super().validate()
        errors = []
        for band1, band2_list in self.bands.items():
            bands = self._split_bands(band2_list)
            if len(set(bands)) != len(bands):
                errors.append(f"value={band2_list} is not a set for key={band1}")
        if errors:
            raise ValueError("\n".join(errors))


class MatchedRefCoaddDiffColorZoomTool(MatchedRefCoaddDiffColorTool):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_color_mmag_zoom_default
        self.produce.metric = NoMetric


class MatchedRefCoaddChiColorTool(MatchedRefCoaddDiffColorTool):
    def setDefaults(self):
        super().setDefaults()
        self.compute_chi = True
        self.produce.plot.yLims = self.limits_chi_default


class MatchedRefCoaddDiffMagTool(MatchedRefCoaddDiffPlot):
    """Tool for diffs between reference and measured coadd mags."""

    mag_y = pexConfig.Field[str](
        default="cmodel_err",
        doc="Flux (magnitude) pexConfig.Field to difference against the x-axis values",
    )
    measure_y_minus_x = pexConfig.Field[bool](
        default=True, doc="Whether to plot the y-axis magnitude minus the x-axis; otherwise x-y if False."
    )

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            # Ensure mag_y is set before any plot finalizes
            self._set_flux_default("mag_y")
            super().finalize()
            self._set_actions()
            name_short_x = self.config_mag_x.name_flux_short
            name_short_y = self.config_mag_y.name_flux_short

            prefix_action = "flux" if self.compute_chi else "mag"
            actionA, actionB = (
                getattr(self.process.buildActions, f"{prefix_action}_{mag}")
                for mag in ((self.mag_y, self.mag_x) if self.measure_y_minus_x else (self.mag_x, self.mag_y))
            )
            action_diff = SubtractVector(actionA=actionA, actionB=actionB)

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
                label_x, label_y = (mag.name_flux for mag in (self.config_mag_x, self.config_mag_y))
                label = f"{label_y} - {label_x}" if self.measure_y_minus_x else f"{label_x} - {label_y}"
                self.produce.plot.yAxisLabel = (
                    f"chi = ({label})/error" if self.compute_chi else f"{label} (mmag)"
                )
            if self.unit is None:
                self.unit = "" if self.compute_chi else "mmag"
            if self.name_prefix is None:
                subtype = "chi" if self.compute_chi else "diff"
                self.name_prefix = f"photom_{name_short_y}_vs_{name_short_x}_mag_{subtype}_{{name_type}}_"
            if self.do_metrics() and not self.produce.metric.units:
                self.produce.metric.units = self.configureMetrics()

    def get_key_flux_y(self) -> str:
        return self.mag_y

    def reconfigure_dependent_magnitudes(
        self,
        key_flux_meas: str | None = None,
        bands_color: dict[str, str] | None = None,
    ):
        if key_flux_meas is not None:
            self.mag_y = key_flux_meas

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_mag_mmag_default


class MatchedRefCoaddDiffMagZoomTool(MatchedRefCoaddDiffMagTool):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_mag_mmag_zoom_default
        self.produce.metric = NoMetric


class MatchedRefCoaddChiMagTool(MatchedRefCoaddDiffMagTool):
    def setDefaults(self):
        super().setDefaults()
        self.compute_chi = True
        self.produce.plot.yLims = self.limits_chi_default


class MatchedRefCoaddDiffPositionTool(MatchedRefCoaddDiffPlot):
    """Tool for diffs between reference and measured coadd positions."""

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
        doc="The key for reference values of the astrometric variable",
        optional=False,
    )
    coord_ref_cos = Field[str](
        doc="The key for reference values of the cosine correction astrometric variable"
        " (i.e. dec if coord_meas is RA)",
        default=None,
        optional=True,
    )
    coord_ref_cos_unit = Field[str](
        doc="astropy unit of coord_ref_cos",
        default="deg",
        optional=True,
    )
    mag_sn = Field[str](default="cmodel_err", doc="Flux (magnitude) field to use for S/N binning on plot")
    # Default coords are in degrees and we want mas
    scale_factor = Field[float](
        doc="The factor to multiply distances by (e.g. the pixel scale if coordinates have pixel units or "
        "the desired spherical coordinate unit conversion otherwise)",
        default=3600000,
    )

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            # Set before MagnitudeScatterPlot.finalize to undo PSF default.
            # Matched ref tables may not have PSF fluxes, or prefer CModel.
            self._set_flux_default("mag_sn")
            super().finalize()
            self._set_actions()
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
                factor = ConstantValue(value=self.scale_factor)
                if self.coord_ref_cos:
                    factor_cos = u.Unit(self.coord_ref_cos_unit).to(u.rad)
                    factor = MultiplyVector(
                        actionA=factor,
                        actionB=CosVector(
                            actionA=MultiplyVector(
                                actionA=ConstantValue(value=factor_cos),
                                actionB=LoadVector(vectorKey=self.coord_meas),
                            )
                        ),
                    )
                self.process.buildActions.diff = MultiplyVector(
                    actionA=factor,
                    actionB=SubtractVector(
                        actionA=self.process.buildActions.pos_meas,
                        actionB=self.process.buildActions.pos_ref,
                    ),
                )
            if self.unit is None:
                self.unit = "" if self.compute_chi else "mas"
            if self.name_prefix is None:
                subtype = "chi" if self.compute_chi else "diff"
                coord_prefix = "" if "coord" in self.coord_meas else "coord_"
                self.name_prefix = (
                    f"astrom_{name_short_y}_vs_{name_short_x}_{coord_prefix}{self.coord_meas}_{subtype}"
                    f"_{{name_type}}_"
                )
            if self.do_metrics() and not self.produce.metric.units:
                self.produce.metric.units = self.configureMetrics()
            if not self.produce.plot.yAxisLabel:
                label = f"({name_short_y} - {name_short_x})"
                coord_suffix = "" if "coord" in name else " coord"
                self.produce.plot.yAxisLabel = (
                    f"chi = ({label} {name}{coord_suffix})/error"
                    if self.compute_chi
                    else f"{label} {name}{coord_suffix} ({self.unit})"
                )

    def get_key_flux_y(self) -> str:
        return self.mag_sn

    def reconfigure_dependent_magnitudes(
        self,
        key_flux_meas: str | None = None,
        bands_color: dict[str, str] | None = None,
    ):
        if key_flux_meas is not None:
            self.mag_sn = key_flux_meas

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_pos_mas_default


class MatchedRefCoaddDiffPositionZoomTool(MatchedRefCoaddDiffPositionTool):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_pos_mas_zoom_default
        self.produce.metric = NoMetric


class MatchedRefCoaddDiffCoordRaTool(MatchedRefCoaddDiffPositionTool):
    def setDefaults(self):
        super().setDefaults()
        self.coord_meas = "coord_ra"
        self.coord_ref = "ref_ra"
        self.coord_ref_cos = "ref_dec"


class MatchedRefCoaddDiffCoordRaZoomTool(MatchedRefCoaddDiffCoordRaTool):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_pos_mas_zoom_default
        self.produce.metric = NoMetric


class MatchedRefCoaddChiCoordRaTool(MatchedRefCoaddDiffCoordRaTool):
    def setDefaults(self):
        super().setDefaults()
        self.compute_chi = True
        self.produce.plot.yLims = self.limits_chi_default


class MatchedRefCoaddDiffCoordDecTool(MatchedRefCoaddDiffPositionTool):
    def setDefaults(self):
        super().setDefaults()
        self.coord_meas = "coord_dec"
        self.coord_ref = "ref_dec"


class MatchedRefCoaddDiffCoordDecZoomTool(MatchedRefCoaddDiffCoordDecTool):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = self.limits_diff_pos_mas_zoom_default
        self.produce.metric = NoMetric


class MatchedRefCoaddChiCoordDecTool(MatchedRefCoaddDiffCoordDecTool):
    def setDefaults(self):
        super().setDefaults()
        self.compute_chi = True
        self.produce.plot.yLims = self.limits_chi_default


class MatchedRefCoaddDiffDistanceTool(MatchedRefCoaddDiffPlot):
    """Tool for distances between matched reference and measured coadd
    objects."""

    key_dist = Field[str](default="match_distance", doc="Distance field key")
    key_dist_err = Field[str](default="match_distanceErr", doc="Distance error field key")
    mag_sn = Field[str](default="cmodel_err", doc="Flux (magnitude) field to use for S/N binning on plot")
    # Default coords are in degrees and we want mas
    scale_factor = Field[float](
        doc="The factor to multiply distances by (e.g. the pixel scale if coordinates have pixel units or "
        "the desired spherical coordinate unit conversion otherwise)",
        default=3600000,
    )

    def finalize(self):
        # Check if it has already been finalized
        if not hasattr(self.process.buildActions, "diff"):
            # Set before MagnitudeScatterPlot.finalize to undo PSF default.
            # Matched ref tables may not have PSF fluxes, or prefer CModel.
            self._set_flux_default("mag_sn")
            super().finalize()
            self._set_actions()

            name_short_x = self.config_mag_x.name_flux_short
            name_short_y = self.config_mag_y.name_flux_short

            self.process.buildActions.dist = LoadVector(vectorKey=self.key_dist)
            if self.compute_chi:
                self.process.buildActions.diff = DivideVector(
                    actionA=self.process.buildActions.dist,
                    actionB=LoadVector(vectorKey=self.key_dist_err),
                )
            else:
                self.process.buildActions.diff = MultiplyVector(
                    actionA=ConstantValue(value=self.scale_factor),
                    actionB=self.process.buildActions.dist,
                )
            if self.unit is None:
                self.unit = "" if self.compute_chi else "mas"
            if self.name_prefix is None:
                subtype = "chi" if self.compute_chi else "diff"
                self.name_prefix = f"astrom_dist_{{name_type}}_{subtype}_"
                self.name_prefix = f"astrom_{name_short_y}_vs_{name_short_x}_dist_{subtype}_{{name_type}}_"
            if self.do_metrics() and not self.produce.metric.units:
                self.produce.metric.units = self.configureMetrics()
            if not self.produce.plot.yAxisLabel:
                label = f"({name_short_y} - {name_short_x}) distance"
                self.produce.plot.yAxisLabel = (
                    f"chi = {label}/error" if self.compute_chi else f"{label} ({self.unit})"
                )

    def get_key_flux_y(self) -> str:
        return self.mag_sn

    def reconfigure_dependent_magnitudes(
        self,
        key_flux_meas: str | None = None,
        bands_color: dict[str, str] | None = None,
    ):
        if key_flux_meas is not None:
            self.mag_sn = key_flux_meas

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = [0, self.limits_diff_pos_mas_default[1]]


class MatchedRefCoaddChiDistanceTool(MatchedRefCoaddDiffDistanceTool):
    def setDefaults(self):
        super().setDefaults()
        self.compute_chi = True
        self.produce.plot.yLims = self.limits_chi_default


class MatchedRefCoaddDiffDistanceZoomTool(MatchedRefCoaddDiffDistanceTool):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.yLims = [0, self.limits_diff_pos_mas_zoom_default[1]]
        self.produce.metric = NoMetric


def reconfigure_diff_matched_defaults(
    config: AnalysisBaseConfig | None = None,
    context: str | None = None,
    key_flux_meas: str | None = None,
    bands_color: dict[str, str] | list[str] | None = None,
    use_any: bool | None = None,
    use_galaxies: bool | None = None,
    use_stars: bool | None = None,
):
    """Reconfigure the default values for config fields of MatchedRefCoaddTool
       and all of its subclasses.

    Parameters
    ----------
    config
        An existing analysis config. Overrides will be applied to any of its
        member MatchedRefCoaddTool atools.
    context
        The context to set. Must be a valid choice for
        MatchedRefCoaddTool.context.
    key_flux_meas
        The key of the measured flux config to use, e.g. "psf". If the key is
        not found, it will search for f"{key}_err", the default name for
        configurations that load error keys as well as fluxes.
    bands_color
        A dictionary keyed by band of comma-separated bands to measure
        colors for, where the color is (key - value). If a list is passed,
        tools will modify the defaults to select only those bands within
        the list (which should also be a set).
    use_any
        Whether to compute metrics for objects of all types.
    use_galaxies
        Whether to compute metrics and plot lines for galaxies only.
    use_stars
        Whether to compute metrics and plot lines for stars only.

    Notes
    -----
    Any kwargs set to None will not change the relevant config field defaults.
    """
    if key_flux_meas is not None:
        keys_flux = tuple(MagnitudeTool.fluxes_default.toDict().keys())
        if key_flux_meas not in keys_flux:
            key_flux_err = f"{key_flux_meas}_err"
            if key_flux_err not in keys_flux:
                raise ValueError(
                    f"{key_flux_meas=} and {key_flux_err} not found in available keys: {keys_flux}"
                    f" (from MagnitudeTool.fluxes_default.toDict().keys())"
                )
            key_flux_meas = key_flux_err

    if context != "custom":
        # These are class attributes and don't need to be changed in subclasses
        # These may end up being changed multiple times with repeated calls,
        # but there isn't a good way to avoid that.
        MagnitudeTool.fluxes_default.ref_matched.name_flux = "True"
        MagnitudeTool.fluxes_default.ref_matched.name_flux_short = "true"
        MagnitudeTool.fluxes_default.ref_matched.key_flux = "ref_{band}_flux"

    def all_subclasses(cls):
        return set(cls.__subclasses__()).union([s for c in cls.__subclasses__() for s in all_subclasses(c)])

    subclasses = all_subclasses(MatchedRefCoaddTool)

    """
    Further context on these kwargs (get it?) and why they're not contexts:

    context (DC2, source_injection, etc)
        This could be made an AnalysisContext, but ChoiceField has the
        benefits of automatic validation. Also, subclasses refer to this config
        field without having to implement separate context functions.
    key_flux_meas
        This is the key to a FluxConfig. Default FluxConfigs could be mapped
        onto an AnalysisContext instead.
    bands_color
        This applies only to color tools and is intended to be set by
        obs package config overrides, e.g. to drop u-band colours. There is no
        way for obs packages to change Tool instance values and the bands
        config field is part of the PipelineTask and not accessible to tools,
        so no obvious alternative exists.
    use_*
        Like context, these apply to all subclasses, but are independent
        booleans rather than exclusive choices.
    """

    # This sets defaults for all known subclasses
    for tool in subclasses:
        tool.reconfigure(
            tool,
            context=context,
            key_flux_meas=key_flux_meas,
            bands_color=bands_color,
            use_any=use_any,
            use_galaxies=use_galaxies,
            use_stars=use_stars,
        )

    # This sets defaults for all existing tools
    # If a pipeline A imports a pipeline B, any atools already set in B will
    # be instantiated before overrides from A are applied. Therefore, changing
    # only the defaults will have no effect on those existing tools.
    if config is not None:
        for tool in config.atools:
            if isinstance(tool, MatchedRefCoaddTool):
                tool: MatchedRefCoaddTool = tool
                tool.reconfigure(
                    context=context,
                    key_flux_meas=key_flux_meas,
                    bands_color=bands_color,
                    use_any=use_any,
                    use_galaxies=use_galaxies,
                    use_stars=use_stars,
                )
