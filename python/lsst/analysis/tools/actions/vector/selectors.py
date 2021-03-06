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
    "FlagSelector",
    "CoaddPlotFlagSelector",
    "SnSelector",
    "ExtendednessSelector",
    "SkyObjectSelector",
    "StarSelector",
    "GalaxySelector",
    "UnknownSelector",
    "VectorSelector",
)

from typing import Optional, cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction


class FlagSelector(VectorAction):
    """The base flag selector to use to select valid sources for QA"""

    selectWhenFalse = ListField[str](
        doc="Names of the flag columns to select on when False", optional=False, default=[]
    )

    selectWhenTrue = ListField[str](
        doc="Names of the flag columns to select on when True", optional=False, default=[]
    )

    def getInputSchema(self) -> KeyedDataSchema:
        allCols = list(self.selectWhenFalse) + list(self.selectWhenTrue)
        return ((col, Vector) for col in allCols)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        """Select on the given flags
        Parameters
        ----------
        table : `Tabular`
        Returns
        -------
        result : `Vector`
            A mask of the objects that satisfy the given
            flag cuts.
        Notes
        -----
        Uses the columns in selectWhenFalse and
        selectWhenTrue to decide which columns to
        select on in each circumstance.
        """
        if not self.selectWhenFalse and not self.selectWhenTrue:
            raise RuntimeError("No column keys specified")
        results: Optional[Vector] = None

        for flag in self.selectWhenFalse:  # type: ignore
            temp = np.array(data[flag.format(**kwargs)] == 0)
            if results is not None:
                results &= temp  # type: ignore
            else:
                results = temp

        for flag in self.selectWhenTrue:
            temp = np.array(data[flag.format(**kwargs)] == 1)
            if results is not None:
                results &= temp  # type: ignore
            else:
                results = temp
        # The test at the beginning assures this can never be None
        return cast(Vector, results)


class CoaddPlotFlagSelector(FlagSelector):
    bands = ListField[str](
        doc="The bands to apply the flags in, takes precedence if band supplied in kwargs",
        default=["i"],
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from super().getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result: Optional[Vector] = None
        bands: tuple[str, ...]
        match kwargs:
            case {"band": band}:
                bands = (band,)
            case {"bands": bands} if not self.bands:
                bands = bands
            case _ if self.bands:
                bands = tuple(self.bands)
            case _:
                bands = ("",)
        for band in bands:
            temp = super().__call__(data, **(kwargs | dict(band=band)))
            if result is not None:
                result &= temp  # type: ignore
            else:
                result = temp
        return cast(Vector, result)

    def setDefaults(self):
        self.selectWhenFalse = [
            "{band}_psfFlux_flag",
            "{band}_pixelFlags_saturatedCenter",
            "{band}_extendedness_flag",
            "xy_flag",
        ]
        self.selectWhenTrue = ["detect_isPatchInner", "detect_isDeblendedSource"]


class VisitPlotFlagSelector(FlagSelector):
    """Select on a set of flags appropriate for making visit-level plots
    (i.e., using sourceTable_visit catalogs).
    """

    def getInputSchema(self) -> KeyedDataSchema:
        yield from super().getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result: Optional[Vector] = None
        temp = super().__call__(data, **kwargs)
        if result is not None:
            result &= temp  # type: ignore
        else:
            result = temp

        return result

    def setDefaults(self):
        self.selectWhenFalse = [
            "psfFlux_flag",
            "pixelFlags_saturatedCenter",
            "extendedness_flag",
            "centroid_flag",
        ]


class SnSelector(VectorAction):
    """Selects points that have S/N > threshold in the given flux type"""

    fluxType = Field[str](doc="Flux type to calculate the S/N in.", default="{band}_psfFlux")
    threshold = Field[float](doc="The S/N threshold to remove sources with.", default=500.0)
    maxSN = Field[float](doc="Maximum S/N to include in the sample (to allow S/N ranges).", default=1e6)
    uncertaintySuffix = Field[str](
        doc="Suffix to add to fluxType to specify uncertainty column", default="Err"
    )
    bands = ListField[str](
        doc="The bands to apply the signal to noise cut in." "Takes precedence if bands passed to call",
        default=[],
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield (fluxCol := self.fluxType), Vector
        yield f"{fluxCol}{self.uncertaintySuffix}", Vector

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        """Makes a mask of objects that have S/N greater than
        self.threshold in self.fluxType
        Parameters
        ----------
        df : `Tabular`
        Returns
        -------
        result : `Vector`
            A mask of the objects that satisfy the given
            S/N cut.
        """
        mask: Optional[Vector] = None
        bands: tuple[str, ...]
        match kwargs:
            case {"band": band}:
                bands = (band,)
            case {"bands": bands} if not self.bands:
                bands = bands
            case _ if self.bands:
                bands = tuple(self.bands)
            case _:
                bands = ("",)
        for band in bands:
            fluxCol = self.fluxType.format(**(kwargs | dict(band=band)))
            errCol = f"{fluxCol}{self.uncertaintySuffix.format(**kwargs)}"
            vec = cast(Vector, data[fluxCol]) / data[errCol]
            temp = (vec > self.threshold) & (vec < self.maxSN)
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp

        # It should not be possible for mask to be a None now
        return np.array(cast(Vector, mask))


class SkyObjectSelector(FlagSelector):
    """Selects sky objects in the given band(s)"""

    bands = ListField[str](
        doc="The bands to apply the flags in, takes precedence if band supplied in kwargs",
        default=["i"],
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from super().getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result: Optional[Vector] = None
        bands: tuple[str, ...]
        match kwargs:
            case {"band": band}:
                bands = (band,)
            case {"bands": bands} if not self.bands:
                bands = bands
            case _ if self.bands:
                bands = tuple(self.bands)
            case _:
                bands = ("",)
        for band in bands:
            temp = super().__call__(data, **(kwargs | dict(band=band)))
            if result is not None:
                result &= temp  # type: ignore
            else:
                result = temp
        return cast(Vector, result)

    def setDefaults(self):
        self.selectWhenFalse = [
            "{band}_pixelFlags_edge",
        ]
        self.selectWhenTrue = ["sky_object"]


class SkySourceSelector(FlagSelector):
    """Selects sky sources from sourceTables"""

    def getInputSchema(self) -> KeyedDataSchema:
        yield from super().getInputSchema()

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result: Optional[Vector] = None
        temp = super().__call__(data, **(kwargs))
        if result is not None:
            result &= temp  # type: ignore
        else:
            result = temp
        return result

    def setDefaults(self):
        self.selectWhenFalse = [
            "pixelFlags_edge",
        ]
        self.selectWhenTrue = ["sky_source"]


class ExtendednessSelector(VectorAction):
    vectorKey = Field[str](
        doc="Key of the Vector which defines extendedness metric", default="{band}_extendedness"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        key = self.vectorKey.format(**kwargs)
        return cast(Vector, data[key])


class StarSelector(ExtendednessSelector):
    extendedness_maximum = Field[float](
        doc="Maximum extendedness to qualify as unresolved, inclusive.", default=0.5, dtype=float
    )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return np.array(cast(Vector, (extendedness >= 0) & (extendedness < self.extendedness_maximum)))


class GalaxySelector(ExtendednessSelector):
    extendedness_minimum = Field[float](
        doc="Minimum extendedness to qualify as resolved, not inclusive.", default=0.5
    )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return cast(Vector, extendedness > self.extendedness_minimum)  # type: ignore


class UnknownSelector(ExtendednessSelector):
    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return extendedness == 9


class VectorSelector(VectorAction):
    """Load a boolean vector from KeyedData and return it for use as a
    selector.
    """

    vectorKey = Field[str](doc="Key corresponding to boolean vector to use as a selection mask")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return cast(Vector, data[self.vectorKey.format(**kwargs)])
