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
    "WPerpPSFMetric",
    "WPerpCModelMetric",
    "XPerpPSFMetric",
    "XPerpCModelMetric",
    "YPerpPSFMetric",
    "YPerpCModelMetric",
    "ValidFracColumnMetric",
)

from lsst.pex.config import Field

from ..actions.scalar import FracInRange, FracNan
from ..actions.vector import LoadVector
from ..analysisParts.stellarLocus import WPerpCModel, WPerpPSF, XPerpCModel, XPerpPSF, YPerpCModel, YPerpPSF
from ..interfaces import AnalysisMetric


class WPerpPSFMetric(WPerpPSF, AnalysisMetric):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "wPerp_psfFlux_sigmaMAD": "mmag",
            "wPerp_psfFlux_median": "mmag",
        }


class WPerpCModelMetric(WPerpCModel, AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "wPerp_cmodelFlux_sigmaMAD": "mmag",
            "wPerp_cmodelFlux_median": "mmag",
        }


class XPerpPSFMetric(XPerpPSF, AnalysisMetric):
    parameterizedBand: bool = False

    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "xPerp_psfFlux_sigmaMAD": "mmag",
            "xPerp_psfFlux_median": "mmag",
        }


class XPerpCModelMetric(XPerpCModel, AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "xPerp_cmodelFlux_sigmaMAD": "mmag",
            "xPerp_cmodelFlux_median": "mmag",
        }


class YPerpPSFMetric(YPerpPSF, AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "yPerp_psfFlux_sigmaMAD": "mmag",
            "yPerp_psfFlux_median": "mmag",
        }


class YPerpCModelMetric(YPerpCModel, AnalysisMetric):
    def setDefaults(self):
        super().setDefaults()

        self.produce.units = {  # type: ignore
            "yPerp_cmodelFlux_sigmaMAD": "mmag",
            "yPerp_cmodelFlux_median": "mmag",
        }


class ValidFracColumnMetric(AnalysisMetric):
    """Calculate the fraction of values in a column that have valid
    numerical values (i.e., not NaN), and that fall within the specified
    "reasonable" range for the values.
    """

    vectorKey = Field[str](doc="Key of column to validate", default="psfFlux")

    def visitContext(self) -> None:
        self.process.buildActions.loadVector = LoadVector()
        self.process.buildActions.loadVector.vectorKey = f"{self.vectorKey}"
        self._setActions(f"{self.vectorKey}")

    def coaddContext(self) -> None:
        self.process.buildActions.loadVector = LoadVector()
        self.process.buildActions.loadVector.vectorKey = "{band}_" + f"{self.vectorKey}"
        self._setActions("{band}_" + f"{self.vectorKey}")

        # Need to pass a mapping of new names so the default names get the
        # band prepended. Otherwise, each subsequent band's metric will
        # overwrite the current one.
        self.produce.newNames = {  # type: ignore
            "validFracColumn": "{band}_validFracColumn",
            "nanFracColumn": "{band}_nanFracColumn",
        }

    def _setActions(self, name: str) -> None:
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

        self.produce.units = {  # type: ignore
            "validFracColumn": "percent",
            "nanFracColumn": "percent",
        }
