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

"""This is a module which defines all the implementation details for the
`Context` base class.
"""

__all__ = ("ContextMeta", "Context", "ContextType", "ContextApplier")

from functools import partial, update_wrapper
from typing import TYPE_CHECKING, Callable, Iterable, Union, cast, overload

from lsst.pex.config.configurableActions import ConfigurableActionStruct

if TYPE_CHECKING:
    from ..interfaces import AnalysisAction


class GetterStandin:
    def __init__(self, base: Union[Context, ContextMeta]):
        self.base = base

    def __call__(self) -> Iterable[ContextMeta]:
        return self.base._contexts


class ContextGetter:
    r"""Return all the individual `Context\ s` that are part of an overall
    `Context`.

    In the case of a single `Context` subclass, this will be a
    set with one element. If the `Context` is a joint `Context` (created
    from more than one individual `Context`\ s, this will be a set of all
    joined `Context`\ s.

    Returns
    -------
    result : `typing.Iterable` of `ContextMeta`
    """

    def __get__(self, instance, klass) -> Callable[..., Iterable[ContextMeta]]:
        if instance is not None:
            return GetterStandin(instance)
        else:
            return GetterStandin(klass)


class ContextApplier:
    @overload
    def __get__(
        self, instance: AnalysisAction, klass: type[AnalysisAction] | None = None
    ) -> Callable[[ContextType], None]: ...

    @overload
    def __get__(self, instance: None, klass: type[AnalysisAction] | None = None) -> ContextApplier: ...

    def __get__(
        self, instance: AnalysisAction | None, klass: type[AnalysisAction] | None = None
    ) -> Callable[[ContextType], None] | ContextApplier:
        if instance is None:
            return self
        part = cast(Callable[[ContextType], None], partial(self.applyContext, instance))
        part = update_wrapper(part, self.applyContext)
        return part

    def __set__(self, instance: AnalysisAction, context: ContextType) -> None:
        self.applyContext(instance, context)

    @staticmethod
    def applyContext(instance: AnalysisAction, context: ContextType) -> None:
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
        # imported here to avoid circular imports at module scope
        from ..interfaces import AnalysisAction

        for ctx in context.getContexts():
            ctx.apply(instance)
        for field in instance._fields:
            match getattr(instance, field):
                case AnalysisAction() as singleField:
                    singleField.applyContext
                    singleField.applyContext(context)
                # type ignore because MyPy is not seeing Pipe_tasks imports
                # correctly (its not formally typed)
                case ConfigurableActionStruct() as multiField:  # type: ignore
                    subField: AnalysisAction
                    for subField in multiField:
                        subField.applyContext(context)


class ContextMeta(type):
    """Metaclass for `Context`, this handles ensuring singleton behavior for
    each `Context`. It also provides the functionality of joining contexts
    together using the | operator.
    """

    _contexts: set[ContextMeta]

    def __new__(cls, *args, **kwargs):
        result = cast(ContextMeta, super().__new__(cls, *args, *kwargs))
        result._contexts = set()
        if result.__name__ != "Context":
            result._contexts.add(result)
        return result

    def apply(cls, action: AnalysisAction) -> None:
        """Apply this context to a given `AnalysisAction`

        This method checks to see if an `AnalysisAction` is aware of this
        `Context`. If it is it calls the actions context method (the class name
        with the first letter lower case)

        Parameters
        ----------
        action : `AnalysisAction`
            The action to apply the `Context` to.
        """
        name = cls.__name__
        name = f"{name[0].lower()}{name[1:]}"
        if hasattr(action, name):
            getattr(action, name)()

    # ignore the conflict with the super type, because we are doing a
    # join
    def __or__(cls, other: ContextMeta | Context) -> Context:  # type: ignore
        """Join multiple Contexts together into a new `Context` instance.

        Parameters
        ----------
        other : `ContextMeta` or `Context`
           #The other `context` to join together with the current `Context`

        Returns
        -------
        jointContext : `Context`
           #A `Context` that is the join of this `Context` and the other.
        """
        if not isinstance(other, (Context, ContextMeta)):
            raise NotImplementedError()
        ctx = Context()
        ctx._contexts = set()
        ctx._contexts |= cls._contexts
        ctx._contexts |= other._contexts
        return ctx

    @staticmethod
    def _makeStr(obj: Union[Context, ContextMeta]) -> str:
        ctxs = list(obj.getContexts())
        if len(ctxs) == 1:
            return ctxs[0].__name__
        else:
            return "|".join(ctx.__name__ for ctx in ctxs)

    def __str__(cls) -> str:
        return cls._makeStr(cls)

    def __repr__(cls) -> str:
        return str(cls)

    getContexts = ContextGetter()


class Context(metaclass=ContextMeta):
    """A Base Context class.

    Instances of this class are used to hold joins of multiple contexts.

    Subclassing this class creates a new independent context.
    """

    _contexts: set[ContextMeta]

    getContexts = ContextGetter()

    def __str__(self) -> str:
        return type(self)._makeStr(self)

    def __repr__(self) -> str:
        return str(self)

    def __or__(self, other: type[Context] | Context) -> Context:
        """Join multiple Contexts together into a new `Context` instance.

        Parameters
        ----------
        other : `ContextMeta` or `Context`
            The other `context` to join together with the current `Context`

        Returns
        -------
        jointContext : `Context`
            A `Context` that is the join of this `Context` and the other.
        """
        if not isinstance(other, (Context, ContextMeta)):
            raise NotImplementedError()
        ctx = Context()
        ctx._contexts = set()
        ctx._contexts |= self._contexts
        ctx._contexts |= other._contexts
        return ctx


ContextType = Union[Context, type[Context]]
"""A type alias to use in contexts where either a Context type or instance
should be accepted (which is most places)
"""
