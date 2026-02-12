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
    "SelectorBase",
    "FlagSelector",
    "CoaddPlotFlagSelector",
    "RangeSelector",
    "SetSelector",
    "SnSelector",
    "ExtendednessSelector",
    "SkyObjectSelector",
    "SkySourceSelector",
    "GoodDiaSourceSelector",
    "StarSelector",
    "GalaxySelector",
    "UnknownSelector",
    "VectorSelector",
    "FiniteSelector",
    "VisitPlotFlagSelector",
    "ThresholdSelector",
    "BandSelector",
    "MatchingFlagSelector",
    "MagSelector",
    "InjectedClassSelector",
    "InjectedGalaxySelector",
    "InjectedObjectSelector",
    "InjectedStarSelector",
    "MatchedObjectSelector",
    "ReferenceGalaxySelector",
    "ReferenceObjectSelector",
    "ReferenceStarSelector",
)

import operator
from typing import Optional, cast

import numpy as np

from lsst.pex.config import Field
from lsst.pex.config.listField import ListField

from ...interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction
from ...math import divide, fluxToMag


class SelectorBase(VectorAction):
    plotLabelKey = Field[str](
        doc="Key to use when populating plot info, ignored if empty string", optional=True, default=""
    )

    def _addValueToPlotInfo(self, value, plotLabelKey=None, **kwargs):
        if "plotInfo" in kwargs:
            if plotLabelKey is not None:
                kwargs["plotInfo"][plotLabelKey] = value
            elif self.plotLabelKey:
                kwargs["plotInfo"][self.plotLabelKey] = value
            else:
                raise RuntimeError(f"No plotLabelKey provided for value {value}, so can't add to plotInfo")


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
            "coord_flag_target",
        ]
        self.selectWhenTrue = ["detect_isPatchInner_target", "detect_isDeblendedSource_target"]

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
            "coord_flag",
            "sky_object",
        ]
        self.selectWhenTrue = ["detect_isPatchInner", "detect_isDeblendedSource"]


class MatchingFlagSelector(CoaddPlotFlagSelector):
    """The default flag selector to apply pre matching.
    The sources are cut down to remove duplicates but
    not on quality.
    """

    def setDefaults(self):
        self.selectWhenFalse = []
        self.selectWhenTrue = ["detect_isPrimary"]


class VisitPlotFlagSelector(FlagSelector):
    """Select on a set of flags appropriate for making visit-level plots
    (i.e., using sourceTable_visit catalogs).
    """

    catalogSuffix = Field[str](doc="The suffix to apply to all the keys.", default="")

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
            "sky_source",
        ]


class RangeSelector(SelectorBase):
    """Selects rows within a range, inclusive of min/exclusive of max."""

    vectorKey = Field[str](doc="Key to select from data")
    maximum = Field[float](doc="The maximum value (exclusive)", default=np.inf)
    minimum = Field[float](doc="The minimum value (inclusive)", default=np.nextafter(-np.inf, 0.0))

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


class SetSelector(SelectorBase):
    """Selects rows with any number of column values within a given set.

    For example, given a set of patches (1, 2, 3), and a set of columns
    (index_1, index_2), return all rows with either index_1 or index_2
     in the set (1, 2, 3).

    Notes
    -----
    The values are given as floats for flexibility. Integers above
    the floating point limit (2^53 + 1 = 9,007,199,254,740,993 for 64 bits)
    will not compare exactly with their float representations.
    """

    vectorKeys = ListField[str](
        doc="Keys to select from data",
        default=[],
        listCheck=lambda x: (len(x) > 0) & (len(x) == len(set(x))),
    )
    values = ListField[float](
        doc="The set of acceptable values",
        default=[],
        listCheck=lambda x: (len(x) > 0) & (len(x) == len(set(x))),
    )

    def getInputSchema(self) -> KeyedDataSchema:
        yield from ((key, Vector) for key in self.vectorKeys)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        """Return a mask of rows with values in the specified set.

        Parameters
        ----------
        data : `KeyedData`

        Returns
        -------
        result : `Vector`
            A mask of the rows with values in the specified set.
        """
        mask = np.zeros_like(data[self.vectorKeys[0]], dtype=bool)
        for key in self.vectorKeys:
            values = cast(Vector, data[key])
            for compare in self.values:
                mask |= values == compare

        return cast(Vector, mask)


class PatchSelector(SetSelector):
    """Select rows within a set of patches."""

    def setDefaults(self):
        super().setDefaults()
        self.vectorKeys = ["patch"]


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
        bandStr = ",".join(bands)
        for band in bands:
            fluxCol = self.fluxType.format(**(kwargs | dict(band=band)))
            fluxInd = fluxCol.find("lux") + len("lux")
            errCol = (
                f"{fluxCol}"[:fluxInd] + f"{self.uncertaintySuffix.format(**kwargs)}" + f"{fluxCol}"[fluxInd:]
            )
            vec = divide(cast(Vector, data[fluxCol]), cast(Vector, data[errCol]))
            temp = (vec > self.threshold) & (vec < self.maxSN)
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp

        plotLabelStr = "({}) > {:.1f}".format(bandStr, self.threshold)
        if self.maxSN < 1e5:
            plotLabelStr += " & < {:.1f}".format(self.maxSN)

        if self.plotLabelKey == "" or self.plotLabelKey is None:
            self._addValueToPlotInfo(plotLabelStr, plotLabelKey="S/N", **kwargs)
        else:
            self._addValueToPlotInfo(plotLabelStr, **kwargs)

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
        super().setDefaults()
        self.selectWhenFalse = [
            "{band}_pixelFlags_edge",
            "{band}_pixelFlags_nodata",
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
        super().setDefaults()
        self.selectWhenFalse = [
            "pixelFlags_edge",
            "pixelFlags_nodata",
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
        super().setDefaults()
        # These default flag names are correct for AP data products
        self.selectWhenFalse = [
            "pixelFlags_bad",
            "pixelFlags_saturatedCenter",
            "pixelFlags_interpolatedCenter",
            "pixelFlags_edge",
            "pixelFlags_nodata",
        ]


class ExtendednessSelector(SelectorBase):
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


class FiniteSelector(VectorAction):
    """Return a mask of finite values for a vector key"""

    vectorKey = Field[str](doc="Key to make a mask of finite values for.")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return cast(Vector, np.isfinite(data[self.vectorKey.format(**kwargs)]))


class VectorSelector(VectorAction):
    """Load a boolean vector from KeyedData and return it for use as a
    selector.
    """

    vectorKey = Field[str](doc="Key corresponding to boolean vector to use as a selection mask")

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.vectorKey, Vector),)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        return cast(Vector, data[self.vectorKey.format(**kwargs)])


class ThresholdSelector(SelectorBase):
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
            mask = np.isin(data[self.vectorKey], bands).flatten()
        else:
            # No band selection is applied, i.e., select all rows
            mask = np.full(len(data[self.vectorKey]), True)  # type: ignore
        return cast(Vector, mask)


class ParentObjectSelector(FlagSelector):
    """Select only parent objects that are not sky objects."""

    def setDefaults(self):
        # This selects all of the parents
        self.selectWhenFalse = [
            "sky_object",
        ]


class ChildObjectSelector(RangeSelector):
    """Select only children from deblended parents"""

    vectorKey = Field[str](doc="Key to select from data", default="parentSourceId")

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
        mask = values > 0

        return cast(Vector, mask)


class MagSelector(SelectorBase):
    """Selects points that have minMag < mag (AB) < maxMag.

    The magnitude is based on the given fluxType.
    """

    fluxType = Field[str](doc="Flux type to calculate the magnitude in.", default="{band}_psfFlux")
    minMag = Field[float](doc="Minimum mag to include in the sample.", default=-1e6)
    maxMag = Field[float](doc="Maximum mag to include in the sample.", default=1e6)
    fluxUnit = Field[str](doc="Astropy unit of flux vector", default="nJy")
    returnMillimags = Field[bool](doc="Use millimags or not?", default=False)
    bands = ListField[str](
        doc="The band(s) to apply the magnitude cut in. Takes precedence if bands passed to call.",
        default=[],
    )

    def getInputSchema(self) -> KeyedDataSchema:
        fluxCol = self.fluxType
        yield fluxCol, Vector

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        """Make a mask of that satisfies self.minMag < mag < self.maxMag.

        The magnitude is based on the flux in self.fluxType.

        Parameters
        ----------
        data : `KeyedData`
            The data to perform the magnitude selection on.

        Returns
        -------
        result : `Vector`
            A mask of the objects that satisfy the given magnitude cut.
        """
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
        bandStr = ",".join(bands)
        for band in bands:
            fluxCol = self.fluxType.format(**(kwargs | dict(band=band)))
            vec = fluxToMag(
                cast(Vector, data[fluxCol]),
                flux_unit=self.fluxUnit,
                return_millimags=self.returnMillimags,
            )
            temp = (vec > self.minMag) & (vec < self.maxMag)
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp

        plotLabelStr = ""
        if self.maxMag < 100:
            plotLabelStr += "({}) < {:.1f}".format(bandStr, self.maxMag)
        if self.minMag > -100:
            if bandStr in plotLabelStr:
                plotLabelStr += " & < {:.1f}".format(self.minMag)
            else:
                plotLabelStr += "({}) < {:.1f}".format(bandStr, self.minMag)
        if self.plotLabelKey == "" or self.plotLabelKey is None:
            self._addValueToPlotInfo(plotLabelStr, plotLabelKey="Mag", **kwargs)
        else:
            self._addValueToPlotInfo(plotLabelStr, **kwargs)

        # It should not be possible for mask to be a None now
        return np.array(cast(Vector, mask))


class InjectedObjectSelector(SelectorBase):
    """A selector for injected objects."""

    vectorKey = Field[str](doc="Key to select from data", default="ref_injected_isPrimary")

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        key = self.vectorKey.format(**kwargs)
        result = cast(Vector, data[key] == 1)
        return result

    def getInputSchema(self) -> KeyedDataSchema:
        yield self.vectorKey, Vector


class InjectedClassSelector(InjectedObjectSelector):
    """A selector for injected objects of a given class."""

    key_class = Field[str](
        doc="Key for the field indicating the class of the object",
        default="ref_source_type",
    )
    key_injection_flag = Field[str](
        doc="Key for the field indicating that the object was not injected" " (per band)",
        default="ref_{band}_injection_flag",
    )
    name_class = Field[str](
        doc="Name of the class of objects",
    )
    value_compare = Field[str](
        doc="Value of the type_key field for objects that are stars",
        default="DeltaFunction",
    )
    value_is_equal = Field[bool](
        doc="Whether the value must equal value_compare to be of this class",
        default=True,
    )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result = super().__call__(data, **kwargs)
        if self.key_injection_flag:
            result &= data[self.key_injection_flag.format(band=kwargs["band"])] == False  # noqa: E712
        values = data[self.key_class]
        result &= (values == self.value_compare) if self.value_is_equal else (values != self.value_compare)
        if self.plotLabelKey:
            self._addValueToPlotInfo(f"injected {self.name_class}", **kwargs)
        return result

    def getInputSchema(self) -> KeyedDataSchema:
        yield from super().getInputSchema()
        yield self.key_class, Vector
        if self.key_injection_flag:
            yield self.key_injection_flag, Vector


class InjectedGalaxySelector(InjectedClassSelector):
    """A selector for injected galaxies."""

    def setDefaults(self):
        self.name_class = "galaxy"
        # Assumes not star == galaxy - if there are injected AGN or other
        # object classes, this will need to be updated
        self.value_is_equal = False


class InjectedStarSelector(InjectedClassSelector):
    """A selector for injected stars."""

    def setDefaults(self):
        self.name_class = "star"


class MatchedObjectSelector(RangeSelector):
    """A selector that selects matched objects with finite distances."""

    def setDefaults(self):
        super().setDefaults()
        self.minimum = 0
        self.vectorKey = "match_distance"


class ReferenceGalaxySelector(ThresholdSelector):
    """A selector that selects galaxies from a catalog with a
    boolean column identifying unresolved sources.
    """

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result = super().__call__(data=data, **kwargs)
        if self.plotLabelKey:
            self._addValueToPlotInfo("reference galaxies", **kwargs)
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

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result = super().__call__(data=data, **kwargs)
        if self.plotLabelKey:
            self._addValueToPlotInfo("reference objects", **kwargs)
        return result

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
            self._addValueToPlotInfo("reference stars", **kwargs)
        return result

    def setDefaults(self):
        super().setDefaults()
        self.op = "eq"
        self.plotLabelKey = "Selection: Stars"
        self.threshold = 1
        self.vectorKey = "refcat_is_pointsource"
