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

__all__ = ("ValidFracColumnMetric",)

from lsst.pex.config import Field

from ..actions.scalar import FracInRange, FracNan
from ..actions.vector import LoadVector
from ..interfaces import AnalysisTool
from .coaddVisit import CoaddVisitConfig


class ValidFracColumnMetric(AnalysisTool, CoaddVisitConfig):
    """Calculate the fraction of values in a column that have valid
    numerical values (i.e., not NaN), and that fall within the specified
    "reasonable" range for the values.
    """

    vectorKey = Field[str](doc="Key of column to validate", default="psfFlux")

    def _setActions(self) -> None:
        if self.context == "coadd":
            name = "{band}_" + f"{self.vectorKey}"
            self.process.buildActions.loadVector = LoadVector()
            self.process.buildActions.loadVector.vectorKey = "{band}_" + f"{self.vectorKey}"

            # Need to pass a mapping of new names so the default names get the
            # band prepended. Otherwise, each subsequent band's metric will
            # overwrite the current one.
            self.produce.metric.newNames = {
                "validFracColumn": "{band}_validFracColumn",
                "nanFracColumn": "{band}_nanFracColumn",
            }
        elif self.context == "visit":
            name = f"{self.vectorKey}"
            self.process.buildActions.loadVector = LoadVector()
            self.process.buildActions.loadVector.vectorKey = f"{self.vectorKey}"
        else:
            raise ValueError(f"Unsupported {self.context=}")

        # The default flux limits of 1e-1 < flux < 1e9 correspond to a
        # magnitude range of 34 < mag < 9 (i.e., reasonably well-measured)
        # objects should all be within this range).
        self.process.calculateActions.validFracColumn = FracInRange(
            vectorKey=name,
            minimum=1.0e-1,
            maximum=1.0e9,
            percent=True,
        )
        self.process.calculateActions.nanFracColumn = FracNan(
            vectorKey=name,
            percent=True,
        )

    def setDefaults(self):
        super().setDefaults()

        self.produce.metric.units = {
            "validFracColumn": "percent",
            "nanFracColumn": "percent",
        }

    def finalize(self):
        AnalysisTool(self).finalize()
        CoaddVisitConfig(self).finalize()
        if not hasattr(self.process.buildActions, "loadVector"):
            self._setActions()
