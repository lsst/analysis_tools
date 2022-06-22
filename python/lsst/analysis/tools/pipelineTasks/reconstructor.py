from __future__ import annotations

from typing import Any

from lsst.daf.butler import Butler, DataId

from lsst.pipe.base.connectionTypes import BaseConnection
from lsst.pipe.base.connections import PipelineTaskConnections

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
