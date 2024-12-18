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

__all__ = ("DiffimMetadataMetricTool",)

from deprecated.sphinx import deprecated
from lsst.pex.config import DictField

from ..actions.scalar import ValueAction
from ..interfaces import AnalysisTool


@deprecated(
    reason=("This tool is superceded by TaskMetadataMetricTool, which should be used instead."),
    version="v29.0",
    category=FutureWarning,
)
class DiffimMetadataMetricTool(AnalysisTool):
    """This tool is designed to extract values from diffim task metadata"""

    metrics = DictField[str, str](doc="The metrics to calculate from the task metadata and their units")

    def finalize(self):
        for name, unit in self.metrics.items():
            setattr(self.process.calculateActions, name, ValueAction(vectorKey=name))
        self.produce.metric.units = dict(self.metrics.items())
