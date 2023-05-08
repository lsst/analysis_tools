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

__all__ = ("ExtendednessTool", "FluxConfig", "MagnitudeTool", "MagnitudeXTool", "SizeConfig", "SizeTool")

import copy

from lsst.pex.config import ChoiceField, Config, ConfigDictField, ConfigField, Field
from lsst.pex.config.configurableActions import ConfigurableActionStructField

from ..actions.vector import (
    CalcShapeSize,
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
    StarSelector,
    VisitPlotFlagSelector,
)
from ..interfaces import AnalysisTool, VectorAction


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
        self.selectors.flagSelector = CoaddPlotFlagSelector()
        self.selectors.flagSelector.bands = ["{band}"]

    def visitContext(self) -> None:
        self.parameterizedBand = False
        self.selectors.flagSelector = VisitPlotFlagSelector()

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

    band_format = Field[str](default="{band}_{key}", doc="Format of band-dependent flux keys.")
    key_flux = Field[str](default=None, doc="Flux field to convert to magnitude on x-axis.")
    key_flux_error = Field[str](default="{key_flux}Err", doc="Flux error field.", optional=True)
    name_flux = Field[str](default=None, doc="Name of the flux/magnitude algorithm/model.")

    def key_flux_band(self, band: str):
        return self.band_format.format(band=band, key=self.key_flux)

    def key_flux_error_band(self, band: str):
        return self.band_format.format(band=band, key=self.key_flux_error.format(key_flux=self.key_flux))


class FluxesDefaultConfig(Config):
    bulge_err = ConfigField[FluxConfig](doc="Bulge model magnitude with errors")
    cmodel_err = ConfigField[FluxConfig](doc="CModel total magnitude with errors")
    disk_err = ConfigField[FluxConfig](doc="Disk model magnitude with errors")
    psf_err = ConfigField[FluxConfig](doc="PSF model magnitude with errors")
    ref_matched = ConfigField[FluxConfig](doc="Reference catalog magnitude")


class MagnitudeTool(AnalysisTool):
    """Compute magnitudes from flux columns.

    Any tool that reads in flux columns and converts them to magnitudes can
    derive from this class and use the _add_flux method to set the
    necessary build actions in their own finalize() methods.
    """

    fluxes_default = FluxesDefaultConfig(
        bulge_err=FluxConfig(key_flux="bdFluxB", name_flux="Bulge"),
        cmodel_err=FluxConfig(key_flux="cModelFlux", name_flux="CModel"),
        disk_err=FluxConfig(key_flux="bdFluxD", name_flux="Disk"),
        psf_err=FluxConfig(key_flux="psfFlux", name_flux="PSF"),
        ref_matched=FluxConfig(
            key_flux="refcat_flux", name_flux="Reference", key_flux_error=None, band_format="{key}_{band}"
        ),
    )

    fluxes = ConfigDictField[str, FluxConfig](  # type: ignore
        default={},
        doc="Flux fields to convert to magnitudes",
    )

    def _add_flux(self, name: str, config: FluxConfig):
        key_flux = config.key_flux_band(band="{band}")
        name_flux = f"flux_{name}"
        self._set_action(self.process.buildActions, name_flux, LoadVector, vectorKey=key_flux)
        if config.key_flux_error is not None:
            # Pre-emptively loaded for e.g. future S/N calculations
            key_flux_err = config.key_flux_error_band(band="{band}")
            self._set_action(
                self.process.buildActions, f"flux_err_{name}", LoadVector, vectorKey=key_flux_err
            )
        self._set_action(self.process.buildActions, f"mag_{name}", ConvertFluxToMag, vectorKey=key_flux)

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

    def _set_flux_default(self, attr):
        """Set own config attr to appropriate string flux name."""
        name_mag = getattr(self, attr)
        # Do nothing if already set - may have been called 2+ times
        if name_mag not in self.fluxes:
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
                self._add_flux(name=name_found, config=value)
            else:
                raise RuntimeError(
                    f"flux={name_mag} not defined in self.fluxes={self.fluxes}"
                    f" and no default configuration found"
                )
            if name_mag != name_found:
                # Essentially appends _err to the name if needed
                setattr(self, attr, name_found)

    def finalize(self):
        super().finalize()
        for key, config in self.fluxes.items():
            self._add_flux(name=key, config=config)


class MagnitudeXTool(MagnitudeTool):
    """A Tool metrics/plots with a magnitude as the dependent variable."""

    mag_x = Field[str](default="", doc="Flux (magnitude) field to bin metrics or plot on x-axis")

    @property
    def config_mag_x(self):
        if self.mag_x not in self.fluxes:
            raise KeyError(f"{self.mag_x=} not in {self.fluxes}; was finalize called?")
        # This is a logic error: it shouldn't be called before finalize
        assert self.mag_x in self.fluxes
        return self.fluxes[self.mag_x]

    def finalize(self):
        super().finalize()
        self._set_flux_default("mag_x")
        key_mag = f"mag_{self.mag_x}"
        subsets = (("xAll", "allSelector"), ("xGalaxies", "galaxySelector"), ("xStars", "starSelector"))
        for name, key in subsets:
            self._set_action(
                self.process.filterActions,
                name,
                DownselectVector,
                vectorKey=key_mag,
                selector=VectorSelector(vectorKey=key),
            )


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


class SizeTool(AnalysisTool):
    """Compute various object size definitions in linear or log space."""

    attr_prefix = Field[str](doc="Prefix to prepend to size names as attrs", default="size_", optional=False)
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
        action = CalcShapeSize(
            colXx=config.key_size.format(suffix="xx"),
            colYy=config.key_size.format(suffix="yy"),
            colXy=config.key_size.format(suffix="xy"),
        )
        return action

    def _get_action_trace(self, config):
        action = CalcShapeSize(
            colXx=config.key_size.format(suffix="xx"), colYy=config.key_size.format(suffix="yy")
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
        self.produce.legendLocation = "lower left"

    def finalize(self):
        super().finalize()
        # A lazy check for whether finalize has already been called
        if hasattr(self.process.filterActions, "yAll"):
            return
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
        self.process.filterActions.yAll = DownselectVector(
            vectorKey=attr, selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey=attr, selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey=attr, selector=VectorSelector(vectorKey="starSelector")
        )
