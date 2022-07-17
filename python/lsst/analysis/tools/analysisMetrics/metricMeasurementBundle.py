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

__all__ = ("MetricMeasurementBundle",)

import json
from collections import UserDict

from lsst.verify import Measurement


class MetricMeasurementBundle(UserDict[str, list[Measurement]]):
    """A specialized dict for storing outputs from multiple `AnalysisMetric`
    actions.

    Keys correspond to the identifier of the action that produced the
    corresponding values. Each value is a list of `~lsst.verift.Measurement`
    objects produced by that `AnalysisMetric` action.

    This object also supports the pydandic interface for serializing to and
    from json compatible objects.
    """

    def json(self):
        result = {}
        for key, value in self.items():
            result[key] = [meas.json for meas in value]
        return json.dumps(result)

    @classmethod
    def parse_obj(cls, data):
        inst = cls()
        for key, value in data.items():
            inst[key] = [Measurement.deserialize(**element) for element in value]
        return inst
