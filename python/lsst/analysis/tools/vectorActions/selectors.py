from __future__ import annotations

__all__ = (
    "FlagSelector",
    "CoaddPlotFlagSelector",
    "SnSelector",
    "ExtendednessSelector",
    "StellarSelector",
    "GalacticSelector",
    "UnknownSelector",
)


from typing import Iterable, Optional, cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField

from ..interfaces import KeyedData, KeyedDataSchema, Vector, VectorAction


class FlagSelector(VectorAction):
    """The base flag selector to use to select valid sources for QA"""

    selectWhenFalse = ListField(
        doc="Names of the flag columns to select on when False", dtype=str, optional=False, default=[]
    )

    selectWhenTrue = ListField(
        doc="Names of the flag columns to select on when True", dtype=str, optional=False, default=[]
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        allCols = list(self.selectWhenFalse) + list(self.selectWhenTrue)  # type: ignore
        return ((col.format(**kwargs), Vector) for col in allCols)

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
            temp = cast(Vector, np.array(data[flag.format(**kwargs)] == 0))
            if results is not None:
                results &= temp  # type: ignore
            else:
                results = temp

        for flag in self.selectWhenTrue:  # type: ignore
            temp = cast(Vector, np.array(data[flag.format(**kwargs)] == 1))
            if results is not None:
                results &= temp  # type: ignore
            else:
                results = temp
        # The test at the beginning assures this can never be None
        return cast(Vector, results)


class CoaddPlotFlagSelector(FlagSelector):
    bands = ListField(
        doc="The bands to apply the flags in, takes precedence if bands supplied in kwargs",
        dtype=str,
        default=["g", "r", "i", "z", "y"],
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        bands = self.bands or kwargs.pop("bands")
        if value := kwargs.pop("band"):
            bands = (value,)
        for band in bands:  # type: ignore
            yield from super().getInputSchema(band=band, **kwargs)

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        result: Optional[Vector] = None
        bands = self.bands or kwargs.pop("bands")
        if value := kwargs.pop("band"):
            bands = (value,)
        for band in bands:  # type: ignore
            temp = super().__call__(data, band=band, **kwargs)
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


class SnSelector(VectorAction):
    """Selects points that have S/N > threshold in the given flux type"""

    fluxType = Field(doc="Flux type to calculate the S/N in.", dtype=str, default="{band}_psfFlux")
    threshold = Field(doc="The S/N threshold to remove sources with.", dtype=float, default=500.0)
    uncertaintySuffix = Field(
        doc="Suffix to add to fluxType to specify uncertainty column", dtype=str, default="Err"
    )
    bands = ListField(
        doc="The bands to apply the signal to noise cut in." "Takes precedence if bands passed to call",
        dtype=str,
        default=["i"],
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        bands = self.bands or kwargs.pop("bands")
        if value := kwargs.pop("band"):
            bands = (value,)
        for band in bands:  # type: ignore
            yield (fluxCol := (cast(str, self.fluxType)).format(**kwargs, band=band)), Vector
            yield f"{fluxCol}{cast(str,self.uncertaintySuffix).format(**kwargs)}", Vector

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
        if not (bands := cast(Iterable[str], self.bands or kwargs.get("bands"))):
            bands = ("",)
        if value := kwargs.pop("band"):
            bands = (value,)
        for band in bands:
            fluxCol = cast(str, self.fluxType).format(**kwargs, band=band)
            errCol = f"{fluxCol}{cast(str,self.uncertaintySuffix).format(**kwargs)}"
            temp = (cast(Vector, data[fluxCol]) / data[errCol]) > cast(float, self.threshold)
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp  # type: ignore

        # It should not be possible for mask to be a None now
        return cast(Vector, mask)


class ExtendednessSelector(VectorAction):
    columnKey = Field(
        doc="Key of the Vector which defines extendedness metric", dtype=str, default="{band}_extendedness"
    )

    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        return ((self.columnKey.format(**kwargs), Vector),)  # type: ignore

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        key = self.columnKey.format(**kwargs)  # type: ignore
        return cast(Vector, data[key])


class StellarSelector(ExtendednessSelector):
    extendedness_maximum = Field(
        doc="Maximum extendedness to qualify as unresolved, inclusive.", default=0.5, dtype=float
    )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return cast(Vector, (extendedness >= 0) & (extendedness < self.extendedness_maximum))


class GalacticSelector(ExtendednessSelector):
    extendedness_minimum = Field(
        doc="Minimum extendedness to qualify as resolved, not inclusive.", default=0.5, dtype=float
    )

    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return cast(Vector, (extendedness >= 0) & (extendedness < self.extendedness_minimum))


class UnknownSelector(ExtendednessSelector):
    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        extendedness = super().__call__(data, **kwargs)
        return extendedness == 9
