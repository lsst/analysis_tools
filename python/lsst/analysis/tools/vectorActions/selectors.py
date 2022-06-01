from __future__ import annotations

__all__ = ("FlagSelector", "CoaddPlotFlagSelector", "SnSelector")


from typing import Iterable, Optional, cast

import numpy as np
from lsst.pex.config import Field
from lsst.pex.config.listField import ListField

from ..interfaces import Tabular, Vector, VectorAction


class FlagSelector(VectorAction):
    """The base flag selector to use to select valid sources for QA"""

    selectWhenFalse = ListField(
        doc="Names of the flag columns to select on when False", dtype=str, optional=False, default=[]
    )

    selectWhenTrue = ListField(
        doc="Names of the flag columns to select on when True", dtype=str, optional=False, default=[]
    )

    def getColumns(self, **kwargs):
        allCols = list(self.selectWhenFalse) + list(self.selectWhenTrue)
        return (col.format(**kwargs) for col in allCols)

    def __call__(self, table: Tabular, **kwargs) -> Vector:
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
            temp = cast(Vector, np.array(table[flag.format(**kwargs)] == 0))
            if results is not None:
                results &= temp  # type: ignore
            else:
                results = temp

        for flag in self.selectWhenTrue:  # type: ignore
            temp = cast(Vector, np.array(table[flag.format(**kwargs)] == 1))
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

    def getColumns(self, **kwargs):
        for band in self.bands or kwargs.get("bands"):  # type: ignore
            yield from super().getColumns(band=band, **kwargs)

    def __call__(self, table: Tabular, **kwargs) -> Vector:
        result: Optional[Vector] = None
        for band in self.bands or kwargs.get("bands"):  # type: ignore
            temp = super.__call__(table, band=band, **kwargs)
            if result is not None:
                result &= temp
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

    fluxType = Field(doc="Flux type to calculate the S/N in.", dtype=str, default="psfFlux")
    threshold = Field(doc="The S/N threshold to remove sources with.", dtype=float, default=500.0)
    uncertantySuffix = Field(
        doc="Suffix to add to fluxType to specify uncertainty column", dtype=str, default="Err"
    )
    bands = ListField(
        doc="The bands to apply the signal to noise cut in." "Takes precedence if bands passed to call",
        dtype=str,
        default=["i"],
    )

    def getColumns(self, **kwargs) -> Iterable[str]:
        if self.bands:
            bands = cast(Iterable[str], self.bands)
        else:
            bands = ("",)
        for band in bands:
            yield (fluxCol := cast(str, self.fluxType)).format(**kwargs, band=band)
            yield f"{fluxCol}_{cast(str,self.uncertantySuffix).format(**kwargs)}"

    def __call__(self, table, **kwargs) -> Vector:
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
        for band in bands:
            fluxCol = cast(str, self.fluxType).format(**kwargs, band=band)
            errCol = f"{fluxCol}_{cast(str,self.uncertantySuffix).format(**kwargs)}"
            temp = cast(Vector, (table[fluxCol] / table[errCol]) > cast(float, self.threshold))
            if mask is not None:
                mask &= temp  # type: ignore
            else:
                mask = temp

        # It should not be possible for mask to be a None now
        return cast(Vector, mask)


class ExtendednessSelector(VectorAction):
    columnKey = Field(
        doc="Key of the column which defines extendedness metric", dtype=str, default="{band}_extendedness_"
    )

    def getColumns(self, **kwargs) -> Iterable[str]:
        return (self.columnKey.format(**kwargs),)  # type: ignore

    def __call__(self, table: Tabular, **kwargs) -> Vector:
        key = self.columnKey.format(**kwargs)  # type: ignore
        return cast(Vector, table[key])


class StellarSelector(ExtendednessSelector):
    extendedness_maximum = Field(
        doc="Maximum extendedness to qualify as unresolved, inclusive.", default=0.5, dtype=float
    )

    def __call__(self, table: Tabular, **kwargs) -> Vector:
        extendedness = super().__call__(table, **kwargs)
        return (extendedness >= 0) & (extendedness < self.extendedness_maximum)


class GalacticSelector(ExtendednessSelector):
    extendedness_minimum = Field(
        doc="Minimum extendedness to qualify as resolved, not inclusive.", default=0.5, dtype=float
    )

    def __call__(self, table: Tabular, **kwargs) -> Vector:
        extendedness = super().__call__(table, **kwargs)
        return (extendedness >= 0) & (extendedness < self.extendedness_minimum)


class UnknownSelector(ExtendednessSelector):
    def __call__(self, table: Tabular, **kwargs) -> Vector:
        extendedness = super().__call__(table, **kwargs)
        return extendedness == 9
