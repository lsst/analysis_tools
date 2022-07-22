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
    "AnalysisAction",
    "KeyedDataAction",
    "ScalarAction",
    "MetricAction",
    "PlotAction",
    "Scalar",
    "KeyedData",
    "KeyedDataTypes",
    "KeyedDataSchema",
    "Vector",
    "AnalysisTool",
    "AnalysisMetric",
    "AnalysisPlot",
)

import warnings
from abc import ABCMeta, abstractmethod
from collections import abc
from numbers import Number
from typing import Any, Iterable, Mapping, MutableMapping, Tuple, Type

import numpy as np
from lsst.pipe.tasks.configurableActions import ConfigurableAction, ConfigurableActionField
from lsst.verify import Measurement
from matplotlib.figure import Figure
from numpy.typing import NDArray

from .contexts import ContextApplier


class ScalarMeta(ABCMeta):
    def __instancecheck__(cls: ABCMeta, instance: Any) -> Any:
        return isinstance(instance, tuple(cls.mro()[1:]))


class Scalar(Number, np.number, metaclass=ScalarMeta):
    """This is an interface only class, and is intended to abstract around all
    the various types of numbers used in Python.

    This has been tried many times with various levels of success in python,
    and this is another attempt. However, as this class is only intended as an
    interface, and not something concrete to use it works.

    Users should not directly instantiate from this class, instead they should
    use a built in python number type, or a numpy number.
    """

    def __init__(self) -> None:
        raise NotImplementedError("Scalar is only an interface and should not be directly instantiated")


Vector = NDArray
"""A Vector is an abstraction around the NDArray interface, things that 'quack'
like an NDArray should be considered a Vector.
"""

KeyedData = MutableMapping[str, Vector | Scalar]
"""KeyedData is an interface where either a `Vector` or `Scalar` can be
retrieved using a key which is of str type.
"""

KeyedDataTypes = MutableMapping[str, Type[Vector] | Type[Number] | Type[np.number]]
r"""A mapping of str keys to the Types which are valid in `KeyedData` objects.
This is useful in conjunction with `AnalysisAction`\ 's ``getInputSchema`` and
``getOutputSchema`` methods.
"""

KeyedDataSchema = Iterable[Tuple[str, Type[Vector] | Type[Number] | Type[np.number]]]
r"""An interface that represents a type returned by `AnalysisAction`\ 's
``getInputSchema`` and ``getOutputSchema`` methods.
"""


class AnalysisAction(ConfigurableAction):
    """Base class interface for the various actions used in analysis tools.

    This extends the basic `ConfigurableAction` class to include interfaces for
    defining what an action expects to consume, and what it expects to produce.
    """

    def __init_subclass__(cls, **kwargs):
        if "getInputSchema" not in dir(cls):
            raise NotImplementedError(f"Class {cls} must implement method getInputSchema")

    # This is a descriptor that functions like a function in most contexts
    # and can be treated as such
    applyContext = ContextApplier()
    r"""Apply a `Context` to an `AnalysisAction` recursively.

    Generally this method is called from within an `AnalysisTool` to
    configure all `AnalysisAction`\ s at one time to make sure that they
    all are consistently configured. However, it is permitted to call this
    method if you are aware of the effects, or from within a specific
    execution environment like a python shell or notebook.

    Parameters
    ----------
    context : `Context`
        The specific execution context, this may be a single context or
        a joint context, see `Context` for more info.
    """

    @abstractmethod
    def getInputSchema(self) -> KeyedDataSchema:
        """Return the schema an `AnalysisAction` expects to be present in the
        arguments supplied to the __call__ method.

        Returns
        -------
        result : `KeyedDataSchema`
            The schema this action requires to be present when calling this
            action, keys are unformatted.
        """
        raise NotImplementedError("This is not implemented on the base class")

    def getOutputSchema(self) -> KeyedDataSchema | None:
        """Return the schema an `AnalysisAction` will produce, if the
        ``__call__`` method returns `KeyedData`, otherwise this may return
        None.

        Returns
        -------
        result : `KeyedDataSchema` or None
            The schema this action will produce when returning from call. This
            will be unformatted if any templates are present. Should return
            None if action does not return `KeyedData`.
        """
        return None

    def getFormattedInputSchema(self, **kwargs) -> KeyedDataSchema:
        """Return input schema, with keys formatted with any arguments supplied
        by kwargs passed to this method.

        Returns
        -------
        result : `KeyedDataSchema`
            The schema this action requires to be present when calling this
            action, formatted with any input arguments (e.g. band='i')
        """
        for key, typ in self.getInputSchema():
            yield key.format_map(kwargs), typ

    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        """Add the supplied inputSchema argument to the class such that it will
        be returned along side any other arguments in a call to
        ``getInputSchema``.

        Parameters
        ----------
        inputSchema : `KeyedDataSchema`
            A schema that is to be merged in with any existing schema when a
            call to ``getInputSchema`` is made.
        """
        warnings.warn(
            f"{type(self)} does not implement adding input schemas, call will do nothing, "
            "this may be expected",
            RuntimeWarning,
        )


class KeyedDataAction(AnalysisAction):
    """A `KeyedDataAction` is an `AnalysisAction` that returns `KeyedData` when
    called.
    """

    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        raise NotImplementedError("This is not implemented on the base class")


class VectorAction(AnalysisAction):
    """A `VectorAction` is an `AnalysisAction` that returns a `Vector` when
    called.
    """

    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        raise NotImplementedError("This is not implemented on the base class")


class ScalarAction(AnalysisAction):
    """A `ScalarAction` is an `AnalysisAction` that returns a `Scalar` when
    called.
    """

    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        raise NotImplementedError("This is not implemented on the base class")

    def getMask(self, **kwargs) -> Vector | slice:
        """Extract a mask if one is passed as key word args, otherwise return
        an empty slice object that can still be used in a getitem call.

        Returns
        -------
        result : `Vector` or `slice`
            The mask passed as a keyword, or a slice object that will return
            a complete Vector when used in getitem.
        """
        if (mask := kwargs.get("mask")) is None:
            mask = slice(None)
        return mask


class MetricAction(AnalysisAction):
    """A `MetricAction` is an `AnalysisAction` that returns a `Measurement` or
    a `Mapping` of `str` to `Measurement` when called.
    """

    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Measurement] | Measurement:
        raise NotImplementedError("This is not implemented on the base class")


class PlotAction(AnalysisAction):
    """A `PlotAction` is an `AnalysisAction` that returns a `Figure` or
    a `Mapping` of `str` to `Figure` when called.
    """

    def getOutputNames(self) -> Iterable[str]:
        """Returns a list of names that will be used as keys if this action's
        call method returns a mapping. Otherwise return an empty Iterable

        Returns
        -------
        result : `Iterable` of `str`
            If a `PlotAction` produces more than one plot, this should be the
            keys the action will use in the returned `Mapping`.
        """
        return tuple()

    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        raise NotImplementedError("This is not implemented on the base class")


class AnalysisTool(AnalysisAction):
    r"""A tool which which calculates a single type of analysis on input data,
    though it may return more than one result.

    AnalysisTools should be used though one of its sub-classes, either an
    `AnalysisMetric` or an `AnalysisPlot`.

    Although `AnalysisTool`\ s are considered a single type of analysis, the
    classes themselves can be thought of as a container. `AnalysisTool`\ s
    are aggregations of `AnalysisAction`\ s to form prep, process, and
    produce stages. These stages allow better reuse of individual
    `AnalysisActions` and easier introspection in contexts such as a notebook
    or interprepter.

    An `AnalysisTool` can be thought of an an individual configuration that
    specifies which `AnalysisAction` should run for each stage.

    The stages themselves are also configurable, allowing control over various
    aspects of the individual `AnalysisAction`\ s.
    """
    prep = ConfigurableActionField[KeyedDataAction](doc="Action to run to prepare inputs")
    process = ConfigurableActionField[AnalysisAction](
        doc="Action to process data into intended form",
    )
    produce = ConfigurableActionField[AnalysisAction](doc="Action to perform any finalization steps")

    parameterizedBand: bool = True
    """Specifies if an `AnalysisTool` may parameterize a band within any field
    in any stage, or if the set of bands is already uniquely determined though
    configuration. I.e. can this `AnalysisTool` be automatically looped over to
    produce a result for multiple bands.
    """

    def __call__(
        self, data: KeyedData, **kwargs
    ) -> Mapping[str, Figure] | Figure | Mapping[str, Measurement] | Measurement:
        bands = kwargs.pop("bands", None)
        if not self.parameterizedBand or bands is None:
            return self._call_single(data, **kwargs)
        results: dict[str, Any] = {}
        if self.identity is not None:
            value_key = f"{{band}}_{self.identity}"
        else:
            value_key = "{band}"
        for band in bands:
            kwargs["band"] = band
            match self._call_single(data, **kwargs):
                case abc.Mapping() as mapping:
                    results.update(mapping.items())
                case value:
                    results[value_key.format(band=band)] = value
        return results

    def _call_single(
        self, data: KeyedData, **kwargs
    ) -> Mapping[str, Figure] | Figure | Mapping[str, Measurement] | Measurement:
        self.populatePrepFromProcess()
        prepped: KeyedData = self.prep(data, **kwargs)  # type: ignore
        processed: KeyedData = self.process(prepped, **kwargs)  # type: ignore
        finalized: Mapping[str, Figure] | Figure | Mapping[str, Measurement] | Measurement = self.produce(
            processed, **kwargs
        )  # type: ignore
        return finalized

    def setDefaults(self):
        super().setDefaults()
        # imported here to avoid circular imports
        from .analysisParts.base import BasePrep, BaseProcess

        self.prep = BasePrep()
        self.process = BaseProcess()

    def getInputSchema(self) -> KeyedDataSchema:
        self.populatePrepFromProcess()
        return self.prep.getInputSchema()

    def populatePrepFromProcess(self):
        """Add additional inputs to the prep stage if supported.

        If the configured prep action supports adding to it's input schema,
        attempt to add the required inputs schema from the process stage to the
        prep stage.

        This method will be a no-op if the prep action does not support this
        feature.
        """
        self.prep.addInputSchema(self.process.getInputSchema())


class AnalysisMetric(AnalysisTool):
    """Specialized `AnalysisTool` for computing metrics.

    The produce stage of `AnalysisMetric` has been specialized such that
    it expects to be assigned to a `MetricAction`, and has a default (set in
    setDefaults) to be `BaseMetricAction`.
    """

    produce = ConfigurableActionField[MetricAction](doc="Action which returns a calculated Metric")

    def setDefaults(self):
        super().setDefaults()
        # imported here to avoid circular imports
        from .analysisParts.base import BaseMetricAction

        self.produce = BaseMetricAction


class AnalysisPlot(AnalysisTool):
    """Specialized `AnalysisTool` for producing plots.

    The produce stage of `AnalysisMetric` has been specialized such that
    it expects to be assigned to a `PlotAction`.
    """

    produce = ConfigurableActionField[PlotAction](doc="Action which returns a plot")

    def getOutputNames(self) -> Iterable[str]:
        """Return the names of the plots produced by this action.

        This will either come from the `PlotAction` if it defines a
        ``getOutputNames`` method (likely if it returns a mapping of figures,
        or a default value is used and a single figure is assumed.

        Returns
        -------
        result : `tuple` of `str`
            Names for each plot produced by this action.
        """
        outNames = tuple(self.produce.getOutputNames())
        if outNames:
            return (f"{self.identity or ''}_{name}" for name in outNames)
        else:
            if self.parameterizedBand:
                return (f"{{band}}_{self.identity or ''}",)
            else:
                return (f"{self.identity or ''}",)
