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
    "Tensor",
    "Scalar",
    "ScalarType",
    "KeyedData",
    "KeyedDataTypes",
    "KeyedDataSchema",
    "Vector",
    "PlotTypes",
    "KeyedResults",
)

from abc import ABCMeta
from collections.abc import Iterable, Mapping, MutableMapping
from numbers import Number
from typing import Any, Protocol, runtime_checkable

import numpy as np
from healsparse import HealSparseMap
from matplotlib.figure import Figure
from numpy.typing import NDArray

from lsst.verify import Measurement


@runtime_checkable
class Tensor(Protocol):
    r"""This is an interface only class and is intended to represent data that
    is 2+ dimensions.

    Technically one could use this for scalars or 1D arrays,
    but for those the Scalar or Vector interface should be preferred.

    `Tensor`\ s abstract around the idea of a multidimensional array, and work
    with a variety of backends including Numpy, CuPy, Tensorflow, PyTorch,
    MXNet, TVM, and mpi4py. This intentionally has a minimum interface to
    comply with the industry standard dlpack which ensures each of these
    backend native types will work.

    To ensure that a `Tensor` is in a desired container (e.g. ndarray) one can
    call the corresponding ``from_dlpack`` method. Whenever possible this will
    be a zero copy action. For instance to work with a Tensor named
    ``input_tensor`` as if it were a numpy object, one would do
    ``image = np.from_dlpack(input_tensor)``.
    """

    ndim: int
    shape: tuple[int, ...]
    strides: tuple[int, ...]

    def __dlpack__(self, /, *, stream: int | None = ...) -> Any: ...

    def __dlpack_device__(self) -> tuple[int, int]: ...


class ScalarMeta(ABCMeta):
    def __instancecheck__(cls: ABCMeta, instance: Any) -> Any:
        return isinstance(instance, tuple(cls.mro()[1:]))


class Scalar(Number, np.number, metaclass=ScalarMeta):  # type: ignore
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


ScalarType = type[Scalar]
"""A type alias for the Scalar interface."""

Vector = NDArray
"""A Vector is an abstraction around the NDArray interface, things that 'quack'
like an NDArray should be considered a Vector.
"""

KeyedData = MutableMapping[str, Vector | Scalar | HealSparseMap | Tensor | Mapping]
"""KeyedData is an interface where either a `Vector`, `Scalar`,
`HealSparseMap`, `Tensor`, or `Mapping` can be retrieved using a key which is
of str type.
"""

KeyedDataTypes = MutableMapping[
    str, type[Vector] | ScalarType | type[HealSparseMap] | type[Tensor] | type[Mapping]
]
r"""A mapping of str keys to the Types which are valid in `KeyedData` objects.
This is useful in conjunction with `AnalysisAction`\ 's ``getInputSchema`` and
``getOutputSchema`` methods.
"""

KeyedDataSchema = Iterable[
    tuple[str, type[Vector] | ScalarType | type[HealSparseMap] | type[Tensor] | type[Mapping]]
]
r"""An interface that represents a type returned by `AnalysisAction`\ 's
``getInputSchema`` and ``getOutputSchema`` methods.
"""

PlotTypes = Figure
"""An interface that represents the various plot types analysis tools supports.
"""

KeyedResults = Mapping[str, PlotTypes | Measurement]
"""A mapping of the return types for an analysisTool."""

type MetricResultType = Mapping[str, Measurement] | Measurement
"""A type alias for the return type of a MetricAction."""

type PlotResultType = Mapping[str, PlotTypes] | PlotTypes
"""A type alias for the return type of a PlotAction."""
