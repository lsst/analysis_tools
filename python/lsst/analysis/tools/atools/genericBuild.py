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
    "ExtendednessTool",
    "FluxConfig",
    "MagnitudeTool",
    "MagnitudeXTool",
    "ObjectClassTool",
    "SizeConfig",
    "SizeTool",
)

import copy

from lsst.pex.config import ChoiceField, Config, ConfigDictField, ConfigField, Field
from lsst.pex.config.configurableActions import ConfigurableActionField, ConfigurableActionStructField

from ..actions.vector import (
    CalcMomentSize,
    ConstantValue,
    ConvertFluxToMag,
    DownselectVector,
    LoadVector,
    Log10Vector,
    MultiplyVector,
    VectorSelector,
)
from ..actions.vector.selectors import (
    CoaddPlotFlagSelector,
    GalaxySelector,
    SelectorBase,
    StarSelector,
    ThresholdSelector,
    VisitPlotFlagSelector,
)
from ..interfaces import AnalysisTool, KeyedData, Vector, VectorAction


class ExtendednessTool(AnalysisTool):
    """Select (non-)extended sources in visit/coadd contexts."""

    extendedness = Field[str](
        default="refExtendedness",
        doc="Extendedness field to select sub-samples with",
    )

    parameterizedBand = Field[bool](
        default=True,
        doc="Does this AnalysisTool support band as a name parameter",
    )

    def coaddContext(self) -> None:
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = ["{band}"]

    def visitContext(self) -> None:
        self.parameterizedBand = False
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()

    def setDefaults(self):
        super().setDefaults()
        # Select any finite extendedness (but still exclude NaNs)
        self.process.buildActions.allSelector = StarSelector(
            vectorKey=self.extendedness, extendedness_maximum=1.0
        )
        self.process.buildActions.galaxySelector = GalaxySelector(vectorKey=self.extendedness)
        self.process.buildActions.starSelector = StarSelector(vectorKey=self.extendedness)


class FluxConfig(Config):
    """Configuration for a flux vector to be loaded and potentially plotted."""

    key_flux = Field[str](default=None, doc="Format of the flux field to convert to magnitudes with {band}.")
    key_flux_error = Field[str](default=None, doc="Format of the flux error field.", optional=True)
    name_flux_short = Field[str](
        default=None,
        doc="Short name of the flux/magnitude algorithm/model to use in metric keys",
    )
    name_flux = Field[str](
        default=None,
        doc="Name of the flux/magnitude algorithm/model to use in plot labels.",
    )

    def key_flux_band(self, band: str):
        return self.key_flux.format(band=band)

    def key_flux_error_band(self, band: str):
        return self.key_flux_error.format(band=band)


class FluxesDefaultConfig(Config):
    bulge_err = ConfigField[FluxConfig](doc="Bulge model magnitude with errors")
    cmodel_err = ConfigField[FluxConfig](doc="CModel total magnitude with errors")
    disk_err = ConfigField[FluxConfig](doc="Disk model magnitude with errors")
    psf_err = ConfigField[FluxConfig](doc="PSF model magnitude with errors")
    ref_matched = ConfigField[FluxConfig](doc="Reference catalog magnitude")


class ObjectSelector(ThresholdSelector):
    """A selector that selects primary objects from an object table."""

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result = super().__call__(data=data, **kwargs)
        if self.plotLabelKey:
            self._addValueToPlotInfo("primary objects", **kwargs)
        return result

    def setDefaults(self):
        super().setDefaults()
        self.op = "eq"
        self.threshold = 1
        self.plotLabelKey = ""
        self.vectorKey = "detect_isPrimary"


class ObjectClassTool(AnalysisTool):
    """Config for tools that compute metrics for multiple classes of object.

    Class refers to e.g. the star-galaxy classification (including other
    types such as AGN), whereas type refers to the more generic (un)resolved
    status of the object. Variability may be included in the future.
    """

    _plurals = {
        # This is a bit of a hack to keep metric names consistent
        # They never changed from all to any, whereas scatterPlot did
        "any": "all",
        "galaxy": "galaxies",
        "star": "stars",
    }

    selection_suffix = Field[str](
        doc="Suffix to append to selector names to summarize selection criteria",
        default="",
    )
    selector_all = ConfigurableActionField[SelectorBase](
        doc="The selector for target objects of all types",
        default=ObjectSelector,
    )
    selector_galaxy = ConfigurableActionField[SelectorBase](
        doc="The selector for target galaxies",
        default=GalaxySelector,
    )
    selector_star = ConfigurableActionField[SelectorBase](
        doc="The selector for target stars",
        default=StarSelector,
    )

    type_any = Field[str](doc="The classification for any type", default="all")
    type_galaxies = Field[str](doc="The classification for galaxies", default="resolved")
    type_stars = Field[str](doc="The classification for galaxies", default="unresolved")

    use_any = Field[bool](doc="Whether any (all types) be a used category of type", default=True)
    use_galaxies = Field[bool](doc="Whether galaxies be a used category of type", default=True)
    use_stars = Field[bool](doc="Whether stars should be a used category of types", default=True)

    def get_class_name_plural(self, object_class: str):
        """Return the plural form of a class name."""
        return self._plurals[object_class]

    def get_class_type(self, object_class: str):
        """Return the type of a given class of objects."""
        match object_class:
            case "any":
                return self.type_any
            case "galaxy":
                return self.type_galaxies
            case "star":
                return self.type_stars

    def get_classes(self):
        """Return all of the classes to be used."""
        classes = []
        if self.use_any:
            classes.append("any")
        if self.use_galaxies:
            classes.append("galaxy")
        if self.use_stars:
            classes.append("star")
        return classes

    def get_name_attr_selector(self, object_class: str, selector_suffix: str = None):
        if selector_suffix is None:
            selector_suffix = self.selection_suffix
        return f"selector{selector_suffix}_{self.get_class_name_plural(object_class)}"

    def get_name_attr_values(self, object_class: str, prefix: str = "y"):
        return f"{prefix}{self.get_class_name_plural(object_class).capitalize()}"

    def get_selector(self, object_class: str):
        """Get the selector for a given object class."""
        match object_class:
            case "any":
                return self.selector_all
            case "galaxy":
                return self.selector_galaxy
            case "star":
                return self.selector_star

    def finalize(self):
        for object_class in self.get_classes():
            name_selector = self.get_name_attr_selector(object_class)
            selector = self.get_selector(object_class)
            # This is a build action because selectors in prep are applied
            # with the and operator. We're not using these to filter all rows
            # but to make several parallel selections.
            setattr(self.process.buildActions, name_selector, selector)

    def setDefaults(self):
        super().setDefaults()
        for selector in (self.selector_galaxy, self.selector_star):
            selector.vectorKey = "refExtendedness"


class MagnitudeTool(ObjectClassTool):
    """Compute magnitudes from flux columns.

    This tool is designed to make it easy to configure multiple fluxes that
    are then converted into magnitudes. For example, a plot might show one
    magnitude on the x-axis, a second on the y-axis, and plot statistics as a
    function of signal-to-noise from a third magnitude.

    The fluxes_default attribute contains commonly used configurations for
    object tables. The "_err" suffix that a flux has an associated error
    column that is needed for some calculation; it can and should be omitted
    if the error column is unneeded. Some flux columns in reference catalogs
    may not have an error at all, such as injection catalogs or truth catalogs
    from simulations.

    Notes
    -----
    Any tool that reads in flux columns and converts them to magnitudes can
    derive from this class and use the _add_flux method to set the
    necessary build actions in their own finalize() methods.
    """

    fluxes_default = FluxesDefaultConfig(
        bulge_err=FluxConfig(
            key_flux="{band}_bdFluxB",
            key_flux_error="{band}_bdFluxBErr",
            name_flux="CModel Bulge",
            name_flux_short="bulge_cModel",
        ),
        cmodel_err=FluxConfig(
            key_flux="{band}_cModelFlux",
            key_flux_error="{band}_cModelFluxErr",
            name_flux="CModel",
            name_flux_short="cModel",
        ),
        disk_err=FluxConfig(
            key_flux="{band}_bdFluxD",
            key_flux_error="{band}_bdFluxDErr",
            name_flux="CModel Disk",
            name_flux_short="disk_cModel",
        ),
        psf_err=FluxConfig(
            key_flux="{band}_psfFlux",
            key_flux_error="{band}_psfFluxErr",
            name_flux="PSF",
            name_flux_short="psf",
        ),
        ref_matched=FluxConfig(
            key_flux="refcat_flux_{band}",
            key_flux_error=None,
            name_flux="Reference",
            name_flux_short="ref",
        ),
    )

    fluxes = ConfigDictField[str, FluxConfig](  # type: ignore
        default={},
        doc="Flux fields to convert to magnitudes",
    )

    def _add_flux(self, name: str, config: FluxConfig, band: str | None = None) -> str:
        """Add requisite buildActions for a given flux.

        Parameters
        ----------
        name
            The name of the flux, without "flux_" prefix.
        config
            The configuration for the flux.
        band
            The name of the band. Default "{band}" assumes the this band is
            the parameterized band.

        Returns
        -------
        name
            The name of the flux, suffixed by band if band is not None.
        """
        if band is None:
            band = "{band}"
        else:
            name = f"{name}_{band}"
        key_flux = config.key_flux_band(band=band)
        name_flux = f"flux_{name}"
        self._set_action(self.process.buildActions, name_flux, LoadVector, vectorKey=key_flux)
        if config.key_flux_error is not None:
            # Pre-emptively loaded for e.g. future S/N calculations
            key_flux_err = config.key_flux_error_band(band=band)
            self._set_action(
                self.process.buildActions, f"flux_err_{name}", LoadVector, vectorKey=key_flux_err
            )
        self._set_action(self.process.buildActions, f"mag_{name}", ConvertFluxToMag, vectorKey=key_flux)
        return name

    def _config_mag(self, name_mag: str = "mag_x"):
        attr = getattr(self, name_mag)
        if attr not in self.fluxes:
            raise KeyError(f"self.{name_mag}={attr} not in {self.fluxes}; was finalize called?")
        return self.fluxes[attr]

    def _finalize_mag(self, name_mag: str = "mag_x", prefix: str = "x"):
        self._set_flux_default(name_mag)
        attr = getattr(self, name_mag)
        key_mag = f"mag_{attr}"
        classes = self.get_classes()
        for object_class in classes:
            name_selector = self.get_name_attr_selector(object_class)
            self._set_action(
                self.process.filterActions,
                self.get_name_attr_values(object_class, prefix=prefix),
                DownselectVector,
                vectorKey=key_mag,
                selector=VectorSelector(vectorKey=name_selector),
            )

    def _set_action(self, target: ConfigurableActionStructField, name: str, action, *args, **kwargs):
        """Set an action attribute on a target tool's struct field.

        Parameters
        ----------
        target
            The ConfigurableActionStructField to set an attribute on.
        name
            The name of the attribute to set.
        action
            The action class to set the attribute to.
        args
            Arguments to pass when initialization the action.
        kwargs
            Keyword arguments to pass when initialization the action.
        """
        if hasattr(target, name):
            attr = getattr(target, name)
            # Setting an attr to a different action is a logic error
            assert isinstance(attr, action)
            # Assert that the action's attributes are identical
            for key, value in kwargs.items():
                if value.__class__.__module__ == "__builtin__":
                    assert getattr(attr, key) == value
        else:
            setattr(target, name, action(*args, **kwargs))

    def _set_flux_default(self, attr: str, band: str | None = None, name_mag: str | None = None) -> str:
        """Set own config attr to appropriate string flux name.

        Parameters
        ----------
        attr
            The name of the attribute to set.
        band
            The name of the band to pass to _add_flux.
        name_mag
            The name of the magnitude to configure. If None, self must already
            have an attr, and name_mag is set to the attr's value.
        """
        name_mag_is_none = name_mag is None
        if name_mag_is_none:
            name_mag = getattr(self, attr)
            complete = name_mag in self.fluxes
        else:
            complete = hasattr(self, attr)
        # Do nothing if already set - may have been called 2+ times
        if not complete:
            name_found = None
            drop_err = False
            # Check if the name with errors is a configured default
            if name_mag.endswith("_err"):
                if hasattr(self.fluxes_default, name_mag):
                    name_found = name_mag
            else:
                if hasattr(self.fluxes_default, name_mag):
                    name_found = name_mag
                # Check if a config with errors exists but not without
                elif hasattr(self.fluxes_default, f"{name_mag}_err"):
                    name_found = f"{name_mag}_err"
                    # Don't load the errors - no _err suffix == unneeded
                    drop_err = True
            if name_found:
                # Copy the config - we don't want to edit in place
                # Other instances may use them
                value = copy.copy(getattr(self.fluxes_default, name_found))
                # Ensure no unneeded error columns are loaded
                if drop_err:
                    value.key_flux_error = None
                self.fluxes[name_found] = value
                name_found = self._add_flux(name=name_found, config=value, band=band)
            else:
                raise RuntimeError(
                    f"flux={name_mag} not defined in self.fluxes={self.fluxes}"
                    f" and no default configuration found"
                )
            if name_mag_is_none and (name_mag != name_found):
                # Essentially appends _err to the name if needed
                setattr(self, attr, name_found)

    def finalize(self):
        super().finalize()
        for key, config in self.fluxes.items():
            self._add_flux(name=key, config=config)


class MagnitudeXTool(MagnitudeTool):
    """A Tool for metrics/plots with a magnitude as the dependent variable."""

    mag_x = Field[str](
        doc="Flux (magnitude) FluxConfig key (in self.fluxes) to bin metrics by or plot on x-axis",
    )

    @property
    def config_mag_x(self):
        return self._config_mag()

    def finalize(self):
        super().finalize()
        self._finalize_mag()


class SizeConfig(Config):
    """Configuration for size vector(s) to be loaded and possibly plotted."""

    has_moments = Field[bool](doc="Whether this size measure is stored as 2D moments.", default=True)
    key_size = Field[str](
        doc="Size column(s) to compute/plot, including moment suffix as {suffix}.",
    )
    log10_size = Field[bool](
        default=True,
        doc="Whether to compute/plot)log10 of the sizes.",
    )
    name_size = Field[str](
        default="size",
        doc="Name of the size (e.g. for axis labels).",
    )
    scale_size = Field[float](
        default=0.2,
        doc="Factor to scale sizes (multiply) by.",
    )
    unit_size = Field[str](
        default="arcsec",
        doc="Unit for sizes.",
    )

    def modify_action(self, action: VectorAction) -> VectorAction:
        if self.log10_size:
            action = Log10Vector(actionA=action)
        return action


class SizeDefaultConfig(Config):
    bulge = ConfigField[SizeConfig](doc="Bulge model size config.")
    disk = ConfigField[SizeConfig](doc="Disk model size config.")
    moments = ConfigField[SizeConfig](doc="Second moments size config.")
    shape_slot = ConfigField[SizeConfig](doc="Shape slot size config.")


class MomentsConfig(Config):
    """Configuration for moment field suffixes."""

    xx = Field[str](doc="Suffix for the x/xx moments.", default="xx")
    xy = Field[str](doc="Suffix for the rho value/xy moments.", default="xy")
    yy = Field[str](doc="Suffix for the y/yy moments.", default="yy")


class SizeTool(ObjectClassTool):
    """Compute various object size definitions in linear or log space.

    This tool is designed to make it easy to configure a size based on the
    definition of the size and the configuration of the columns that it is
    read from.

    It currently only supports a single size but may be refactored to
    support arbitrary sizes, like MagnitudeTool.
    """

    attr_prefix = Field[str](doc="Prefix to prepend to size names as attrs", default="size_", optional=False)
    config_moments = ConfigField[MomentsConfig](
        doc="Configuration for moment field names", default=MomentsConfig
    )
    is_covariance = Field[bool](
        doc="Whether this size has multiple fields as for a covariance matrix."
        " If False, the XX/YY/XY terms are instead assumed to map to sigma_x/sigma_y/rho.",
        default=True,
    )
    sizes_default = SizeDefaultConfig(
        bulge=SizeConfig(key_size="{band}_bdReB", name_size="CModel Bulge $R_{eff}$", has_moments=False),
        disk=SizeConfig(key_size="{band}_bdReD", name_size="CModel Disk $R_{eff}$", has_moments=False),
        moments=SizeConfig(key_size="{band}_i{suffix}", name_size="Second moment radius"),
        shape_slot=SizeConfig(key_size="shape_{suffix}", name_size="Shape slot radius"),
    )
    size_type = ChoiceField[str](
        doc="The type of size to calculate",
        allowed={
            "determinantRadius": "The (matrix) determinant radius from x/y moments.",
            "traceRadius": "The (matrix) trace radius from x/y moments.",
            "singleColumnSize": "A pre-computed size from a single column.",
        },
        optional=False,
    )
    size_y = Field[str](default=None, doc="Name of size field to plot on y axis.")
    sizes = ConfigDictField[str, SizeConfig](  # type: ignore
        default={},
        doc="Size fields to add to build actions",
    )

    def _check_attr(self, name_size: str):
        """Check if a buildAction has already been set."""
        attr = self.get_attr_name(name_size)
        if hasattr(self.process.buildActions, attr):
            raise RuntimeError(f"Can't re-set size build action with already-used {attr=} from {name_size=}")

    def _get_action_determinant(self, config):
        action = CalcMomentSize(
            colXx=config.key_size.format(suffix=self.config_moments.xx),
            colYy=config.key_size.format(suffix=self.config_moments.yy),
            colXy=config.key_size.format(suffix=self.config_moments.xy),
            is_covariance=self.is_covariance,
        )
        return action

    def _get_action_trace(self, config):
        action = CalcMomentSize(
            colXx=config.key_size.format(suffix=self.config_moments.xx),
            colYy=config.key_size.format(suffix=self.config_moments.yy),
            is_covariance=self.is_covariance,
        )
        return action

    def _get_action_single_column(self, config):
        action = LoadVector(vectorKey=config.key_size)
        return action

    def get_attr_name(self, name_size):
        """Return the build action attribute for a size of a given name."""
        return f"{self.attr_prefix}{name_size}"

    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.legendLocation = "lower left"

    def finalize(self):
        # A lazy check for whether finalize has already been called
        classes = self.get_classes()
        if hasattr(self.process.filterActions, self.get_name_attr_values(classes[0])):
            return
        super().finalize()
        if not self.size_y:
            raise ValueError("Must specify size_y")
        elif self.size_y not in self.sizes:
            if size_y := getattr(self.sizes_default, self.size_y, None):
                self.sizes[self.size_y] = size_y
            else:
                raise RuntimeError(f"{self.size_y=} not found in {self.sizes=} or {self.sizes_default=}")

        if self.size_type == "determinantRadius":
            get_action = self._get_action_determinant
        elif self.size_type == "traceRadius":
            get_action = self._get_action_trace
        elif self.size_type == "singleColumnSize":
            get_action = self._get_action_single_column
        else:
            raise ValueError(f"Unsupported {self.size_type=}")

        for name, config in self.sizes.items():
            self._check_attr(name)
            action = config.modify_action(
                MultiplyVector(
                    actionA=get_action(config=config),
                    actionB=ConstantValue(value=config.scale_size),
                )
            )
            setattr(self.process.buildActions, self.get_attr_name(name), action)

        attr = self.get_attr_name(self.size_y)
        classes = self.get_classes()
        for object_class in classes:
            setattr(
                self.process.filterActions,
                self.get_name_attr_values(object_class),
                DownselectVector(
                    vectorKey=attr,
                    selector=VectorSelector(vectorKey=self.get_name_attr_selector(object_class)),
                ),
            )
