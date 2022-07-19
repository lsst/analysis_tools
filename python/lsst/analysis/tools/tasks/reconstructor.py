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

from typing import Any

from lsst.daf.butler import Butler, DataId
from lsst.pipe.base.connections import PipelineTaskConnections
from lsst.pipe.base.connectionTypes import BaseConnection

from .base import AnalysisBaseConfig


def reconstructAnalysisTools(
    butler: Butler, collection: str, label: str, dataId: DataId
) -> tuple[AnalysisBaseConfig, dict[str, Any]]:
    configDSType = f"{label}_config"
    config = butler.get(configDSType, collections=(collection,))

    connections: PipelineTaskConnections = config.connections.ConnectionsClass(config=config)
    inputs = {}

    for name in connections.inputs:
        connection: BaseConnection = getattr(connections, name)
        dsName = connection.name
        inputs[name] = butler.get(dsName, dataId=dataId, collections=(collection,))

    return (config, inputs)
