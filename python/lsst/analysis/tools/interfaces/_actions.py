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
    "VectorAction",
    "ScalarAction",
    "MetricResultType",
    "MetricAction",
    "PlotResultType",
    "PlotAction",
    "JointResults",
    "JointAction",
)

import warnings
from abc import abstractmethod
from dataclasses import dataclass
from typing import Iterable

from lsst.pex.config.configurableActions import ConfigurableAction, ConfigurableActionField

from ..contexts import ContextApplier
from ._interfaces import KeyedData, KeyedDataSchema, MetricResultType, PlotResultType, Scalar, Vector


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
    def __call__(self, data: KeyedData, **kwargs) -> MetricResultType:
        raise NotImplementedError("This is not implemented on the base class")


class PlotAction(AnalysisAction):
    """A `PlotAction` is an `AnalysisAction` that returns a `PlotType` or
    a `Mapping` of `str` to `PlotType` when called.
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
    def __call__(self, data: KeyedData, **kwargs) -> PlotResultType:
        raise NotImplementedError("This is not implemented on the base class")


class NoPlot(PlotAction):
    """This is a sentinel class to indicate that there is no plotting action"""


class NoMetric(MetricAction):
    """This is a sentinel class to indicate that there is no Metric action"""


@dataclass
class JointResults:
    plot: PlotResultType | None
    metric: MetricResultType | None


class JointAction(AnalysisAction):
    """A `JointAction` is an `AnalysisAction` that is a composite of a
    `PlotAction` and a `MetricAction`
    """

    metric = ConfigurableActionField[MetricAction](doc="Action to run that will produce one or more metrics")
    plot = ConfigurableActionField[PlotAction](doc="Action to run that will produce one or more plots")

    def __call__(self, data: KeyedData, **kwargs) -> JointResults:
        if isinstance(self.plot, NoPlot):
            plots = None
        else:
            plots = self.plot(data, **kwargs)
        if isinstance(self.metric, NoMetric):
            metrics = None
        else:
            metrics = self.metric(data, **kwargs)
        return JointResults(plots, metrics)

    def getInputSchema(self) -> KeyedDataSchema:
        yield from self.metric.getInputSchema()
        yield from self.plot.getInputSchema()

    def getOutputNames(self) -> Iterable[str]:
        """Returns a list of names that will be used as keys if this action's
        call method returns a mapping. Otherwise return an empty Iterable

        Returns
        -------
        result : `Iterable` of `str`
            If a `PlotAction` produces more than one plot, this should be the
            keys the action will use in the returned `Mapping`.
        """
        return self.plot.getOutputNames()
