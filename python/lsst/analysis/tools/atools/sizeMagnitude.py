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

import logging

from deprecated.sphinx import deprecated

from ..actions.vector import CoaddPlotFlagSelector, VisitPlotFlagSelector
from .genericBuild import SizeTool
from .genericProduce import MagnitudeScatterPlot

logging.basicConfig()
_LOG = logging.getLogger(__name__)
_LOG.setLevel(logging.WARNING)


class SizeMagnitudePlot(SizeTool, MagnitudeScatterPlot):

    _parameterizedBand: bool = True

    # TODO: Remove the getter and setting in DM-54864.
    # Instead, just convert _parameterizedBand to parameterizedBand.
    @property
    def parameterizedBand(self) -> bool:
        return self._parameterizedBand

    @parameterizedBand.setter
    def parameterizedBand(self, value: bool) -> None:
        _LOG.info(
            "Setting 'parameterizedBand' of a `SizeMagnitudePlot' is a no-op. "
            "If you see this message, it is likely because you are trying "
            "to read older version of configs with newer versions of the "
            "LSST Science Pipelines."
        )

    # TODO: Remove the getter and setter in DM-54864.
    @property
    @deprecated(reason="This is no longer used.", version="v31")
    def extendedness(self) -> None:
        """A config-like attribute for backward compatibility.

        This does not do anything but enable reading old configs.
        """

    @extendedness.setter
    def extendedness(self, value: str) -> None:
        _LOG.info(
            "Setting 'extendedness' of a `SizeMagnitudePlot' is a no-op. "
            "If you see this message, it is likely because you are trying "
            "to read older version of configs with newer versions of the "
            "LSST Science Pipelines."
        )

    def coaddContext(self) -> None:
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = []

    def visitContext(self) -> None:
        self.prep.selectors.flagSelector = VisitPlotFlagSelector()

    def finalize(self):
        # TODO: Investigate why MagnitudeScatterPlot.finalize(self) is called
        # always, even if super().finalize() is omitted
        super().finalize()
        if not self.produce.plot.yAxisLabel:
            size = self.sizes[self.size_y]
            self.produce.plot.yAxisLabel = (
                f"log10({size.name_size}/{size.unit_size})"
                if size.log10_size
                else f"{size.name_size} ({size.unit_size})"
            )
