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

__all__ = ("CountPatches",)

from ..actions.scalar.scalarActions import CountUniqueAction
from ..actions.vector import FiniteSelector, LoadVector
from ..interfaces import AnalysisTool


class CountPatches(AnalysisTool):
    """An atool to count the patches in a tract.
    Counts the patches that make it into the
    objectTable_tract.
    """

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.patchSelector = FiniteSelector(vectorKey="{band}_ra")
        self.process.buildActions.patch = LoadVector(vectorKey="patch")

        self.process.calculateActions.patchCount = CountUniqueAction()
        self.process.calculateActions.patchCount.vectorKey = "patch"

        self.produce.metric.units = {"patchCount": ""}

        self.produce.metric.newNames = {
            "patchCount": "{band}_patchCount",
        }
