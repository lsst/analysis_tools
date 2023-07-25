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
    "RangeSelector",
    "SnSelector",
    "ExtendednessSelector",
    "SkyObjectSelector",
    "SkySourceSelector",
    "GoodDiaSourceSelector",
    "StarSelector",
    "GalaxySelector",
    "UnknownSelector",
    "VectorSelector",
    "VisitPlotFlagSelector",
    "ThresholdSelector",
    "BandSelector",
)

import operator
from typing import Optional, cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction


class SelectorBase(VectorAction):
    plotLabelKey = Field[str](
        doc="Key to use when populating plot info, ignored if empty string", optional=True, default=""
    )

    def _addValueToPlotInfo(self, value, **kwargs):
        if "plotInfo" in kwargs and self.plotLabelKey:
            kwargs["plotInfo"][self.plotLabelKey] = value


class FlagSelector(VectorAction):
    """The base flag selector to use to select valid sources for QA."""

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
        data : `KeyedData`

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
    """This default setting makes it take the band from
    the kwargs.
    """

    bands = ListField[str](
        doc="The bands to apply the flags in, takes precedence if band supplied in kwargs",
        default=[],
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from super().getInputSchema()

    def refMatchContext(self):
        self.selectWhenFalse = [
            "{band}_psfFlux_flag_target",
            "{band}_pixelFlags_saturatedCenter_target",
            "{band}_extendedness_flag_target",
            "xy_flag_target",
        ]
        self.selectWhenTrue = ["detect_isPatchInner_target",
                               "detect_isDeblendedSource_target"]

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result: Optional[Vector] = None
        bands: tuple[str, ...]
        match kwargs:
            case {"band": band} if not self.bands and self.bands == []:
                bands = (band,)
            case {"bands": bands} if not self.bands and self.bands == []:
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
        self.selectWhenTrue = ["detect_isPatchInner",
                               "detect_isDeblendedSource"]


class VisitPlotFlagSelector(FlagSelector):
    """Select on a set of flags appropriate for making visit-level plots
    (i.e., using sourceTable_visit catalogs).
    """

    catalogSuffix = Field[str](
        doc="The suffix to apply to all the keys.", default="")

    def getInputSchema(self) -> KeyedDataSchema:
        yield from super().getInputSchema()

    def refMatchContext(self):
        self.selectWhenFalse = [
            "psfFlux_flag_target",
            "pixelFlags_saturatedCenter_target",
            "extendedness_flag_target",
            "centroid_flag_target",
        ]


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


class RangeSelector(VectorAction):
    """Selects rows within a range, inclusive of min/exclusive of max."""

    vectorKey = Field[str](doc="Key to select from data")
    maximum = Field[float](doc="The maximum value", default=np.Inf)
    minimum = Field[float](doc="The minimum value", default=np.nextafter(-np.Inf, 0.0))

    def getInputSchema(self) -> KeyedDataSchema:
        yield self.vectorKey, Vector

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        """Return a mask of rows with values within the specified range.

        Parameters
        ----------
        data : `KeyedData`

        Returns
        -------
        result : `Vector`
            A mask of the rows with values within the specified range.
        """
        values = cast(Vector, data[self.vectorKey])
        mask = (values >= self.minimum) & (values < self.maximum)

        return cast(Vector, mask)


class SnSelector(SelectorBase):
    """Selects points that have S/N > threshold in the given flux type."""

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
        fluxCol = self.fluxType
        fluxInd = fluxCol.find("lux") + len("lux")
        errCol = f"{fluxCol}"[:fluxInd] + f"{self.uncertaintySuffix}" + f"{fluxCol}"[fluxInd:]
        yield fluxCol, Vector
        yield errCol, Vector

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        """Makes a mask of objects that have S/N greater than
        self.threshold in self.fluxType

        Parameters
        ----------
        data : `KeyedData`
            The data to perform the selection on.

        Returns
        -------
        result : `Vector`
            A mask of the objects that satisfy the given
            S/N cut.
        """

        self._addValueToPlotInfo(self.threshold, **kwargs)
        mask: Optional[Vector] = None
        bands: tuple[str, ...]
        match kwargs:
            case {"band": band} if not self.bands and self.bands == []:
                bands = (band,)
            case {"bands": bands} if not self.bands and self.bands == []:
                bands = bands
            case _ if self.bands:
                bands = tuple(self.bands)
            case _:
                bands = ("",)
        for band in bands:
            fluxCol = self.fluxType.format(**(kwargs | dict(band=band)))
            fluxInd = fluxCol.find("lux") + len("lux")
            errCol = f"{fluxCol}"[:fluxInd] + f"{self.uncertaintySuffix.format(**kwargs)}" + f"{fluxCol}"[fluxInd:]
            vec = cast(Vector, data[fluxCol]) / data[errCol]
            temp = (vec > self.threshold) & (vec < self.maxSN)
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp

        # It should not be possible for mask to be a None now
        return np.array(cast(Vector, mask))


class SkyObjectSelector(FlagSelector):
    """Selects sky objects in the given band(s)."""

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
            case {"band": band} if not self.bands and self.bands == []:
                bands = (band,)
            case {"bands": bands} if not self.bands and self.bands == []:
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
    """Selects sky sources from sourceTables."""

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


class GoodDiaSourceSelector(FlagSelector):
    """Selects good DIA sources from diaSourceTables."""

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
            "base_PixelFlags_flag_bad",
            "base_PixelFlags_flag_suspect",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_interpolated",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_edge",
        ]


class ExtendednessSelector(VectorAction):
    """A selector that picks between extended and point sources."""

    vectorKey = Field[str](
        doc="Key of the Vector which defines extendedness metric", default="{band}_extendedness"
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        key = self.vectorKey.format(**kwargs)
        return cast(Vector, data[key])


class StarSelector(ExtendednessSelector):
    """A selector that picks out stars based off of their
    extendedness values.
    """

    extendedness_maximum = Field[float](
        doc="Maximum extendedness to qualify as unresolved, inclusive.", default=0.5, dtype=float
    )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return np.array(cast(Vector, (extendedness >= 0) & (extendedness < self.extendedness_maximum)))


class GalaxySelector(ExtendednessSelector):
    """A selector that picks out galaxies based off of their
    extendedness values.
    """

    extendedness_minimum = Field[float](
        doc="Minimum extendedness to qualify as resolved, not inclusive.", default=0.5
    )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return cast(Vector, extendedness > self.extendedness_minimum)  # type: ignore


class UnknownSelector(ExtendednessSelector):
    """A selector that picks out unclassified objects based off of their
    extendedness values.
    """

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


class ThresholdSelector(VectorAction):
    """Return a mask corresponding to an applied threshold."""

    op = Field[str](doc="Operator name.")
    threshold = Field[float](doc="Threshold to apply.")
    vectorKey = Field[str](doc="Name of column")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        mask = getattr(operator, self.op)(data[self.vectorKey], self.threshold)
        return cast(Vector, mask)


class BandSelector(VectorAction):
    """Makes a mask for sources observed in a specified set of bands."""

    vectorKey = Field[str](doc="Key of the Vector which defines the band", default="band")
    bands = ListField[str](
        doc="The bands to select. `None` indicates no band selection applied.",
        default=[],
    )

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        bands: Optional[tuple[str, ...]]
        match kwargs:
            case {"band": band} if not self.bands and self.bands == []:
                bands = (band,)
            case {"bands": bands} if not self.bands and self.bands == []:
                bands = bands
            case _ if self.bands:
                bands = tuple(self.bands)
            case _:
                bands = None
        if bands:
            mask = np.in1d(data[self.vectorKey], bands)
        else:
            # No band selection is applied, i.e., select all rows
            mask = np.full(len(data[self.vectorKey]), True)  # type: ignore
        return cast(Vector, mask)
