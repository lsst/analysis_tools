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

from ..analysisParts.stellarLocus import WPerpCModel, WPerpPSF, XPerpCModel, XPerpPSF, YPerpCModel, YPerpPSF
from ..interfaces import AnalysisMetric
from ..actions.scalar import FracInRange
from ..actions.vector import LoadVector, FlagSelector
from ..actions.keyedData import KeyedDataSelectorAction


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

    # parameterizedBand: bool = True

    columnKey: str = "psfFlux"

    def visitContext(self) -> None:
        self.prep.selectors.keyedDataSelector = KeyedDataSelectorAction()
        self._setActions(f"{self.columnKey}")

    def coaddContext(self) -> None:
        self.prep.selectors.keyedDataSelector = KeyedDataSelectorAction()
        # self.process.buildActions.loadVector = LoadVector()
        # self.process.buildActions.loadVector.vectorKey = f"{{band}}_{self.columnKey}"
        self._setActions(f"{{band}}_{self.columnKey}")

        # Need to pass a mapping of new names so the default names get the
        # band prepended. Otherwise, each subsequent band's metric will
        # overwrite the current one.
        self.produce.newNames = {
            "validFracColumn": "{band}_validFracColumn",
        }

    def _setActions(self, name: str) -> None:
        self.process.calculateActions.validFracColumn = FracInRange(
            vectorKey=name,
            minimum=1.0e-2,
            maximum=1.0e6,
            percent=True,
        )

    def setDefaults(self):
        super().setDefaults()

        # self.process.buildActions.rangeSelector = RangeSelector()
        #self.process.buildActions.fracInRangeSelector = FracInRange()

        # select values within range
        #self.process.buildActions.fracInRangeSelector.vectorKey = "psfFlux"
        #self.process.buildActions.fracInRangeSelector.minimum = 1e-2
        #self.process.buildActions.fracInRangeSelector.maximum = 1e6

#        self.process.calculateActions.validFracColumn = FracInRange(
#            vectorKey="psfFlux",
#            minimum=1.0e-2,
#            maximum=1.0e6,
#            percent=True,
#        )

        self.produce.units = {"validFracColumn": "percent"}

        # the final name in the qualification is used as a key to insert
        # the calculation into KeyedData
        #self.process.filterActions.allInRange = DownselectVector(
        #    vectorKey=self.process.buildActions.rangeSelector.vectorKey,
        #    selector=self.process.buildActions.rangeSelector
        #)

        #self.process.calculateActions.ValidFracColumnMetric = CountAction(vectorKey="allInRange")

        #self.produce.units = {"ValidFracColumnMetric": "ct"}

        #self.produce.units = {  # type: ignore
        #    "yPerp_psfFlux_sigmaMAD": "mmag",
        #    "yPerp_psfFlux_median": "mmag",
        #}
