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

__all__ = ("reconstructAnalysisTools", "getPlotDatasetTypeNames")

from typing import TYPE_CHECKING, Any, Callable, Iterable

from lsst.pipe.base.connections import PipelineTaskConnections, iterConnections
from lsst.pipe.base.connectionTypes import BaseConnection

from .base import AnalysisBaseConfig

if TYPE_CHECKING:
    from lsst.daf.butler import Butler, DataId


def reconstructAnalysisTools(
    butler: Butler,
    collection: str,
    label: str,
    dataId: DataId,
    callback: Callable[[dict[str, Any], DataId], dict[str, Any]] | None,
) -> tuple[AnalysisBaseConfig, dict[str, Any]]:
    """Reconstructs the analysis tools used to produce metrics and plots in a
    task and all input data required.

    Parameters
    ----------
    butler : `~lsst.daf.butler.Butler`
        The butler where the data is stored.
    collection : `str`
        Collection within the butler associated with desired data.
    label : `str`
        The label from the `~lsst.pipe.base.Pipeline` associated with the task
        whose tools are to be reconstructed.
    dataId : `~lsst.daf.butler.DataId`
        Identifier for which data to retrieve.
    callback : `~typing.Callable` or None
        An optional function which can transform the data after it has been
        loaded from the butler. The function must take a dict of strings to
        data products, and the DataId. The function must return a dict of
        string to data products. The returned dict is what will be returned
        by `reconstructAnalysisTools`

    Returns
    -------
    config : `AnalysisBaseConfig`
        The configuration of the task run to produce metrics and plots. This
        config contains all the `AnalysisTools` as configured when the task
        produced the data.
    data : `dict` of `str` to `Any`
        The data that went into producing metrics and plots.
    """
    configDSType = f"{label}_config"
    config = butler.get(configDSType, collections=(collection,))

    connections: PipelineTaskConnections = config.connections.ConnectionsClass(config=config)
    inputs: dict[str, Any] = {}

    for name in connections.inputs:
        connection: BaseConnection = getattr(connections, name)
        dsName = connection.name
        # If this is a multiple connection, query the butler for all the
        # inputs for this dataset type name
        if connection.multiple:
            container = []
            for ref in set(
                butler.registry.queryDatasets(
                    dsName, dataId=dataId, findFirst=True, collections=(collection,)
                )
            ):
                container.append(butler.get(ref, collections=(collection,)))
            inputs[name] = container
        else:
            inputs[name] = butler.get(dsName, dataId=dataId, collections=(collection,))

    if callback is not None:
        inputs = callback(inputs, dataId)

    return (config, inputs)


def getPlotDatasetTypeNames(
    butler: Butler,
    collections: str | Iterable[str],
    label: str | None = None,
) -> Iterable[str]:
    """Get the dataset type names for plots (anything with StorageClass="Plot")
    from butler collections.

    Parameters
    ----------
    butler : `~lsst.daf.butler.Butler`
        The butler where the data is stored.
    collections : `str` or `list` [`str`]
        Collections within the butler to query for datasets containing plots.
    label : `str`, optional
        The label from the `~lsst.pipe.base.Pipeline` associated with the task
        whose plots are to be queried. If no label is given, all requested
        collections will be queried.

    Returns
    -------
    plotNames : `list` [`str`]
        Plot dataset type names.
    """
    if label is not None:
        configs = [butler.get(f"{label}_config", collections=collections)]
    else:
        configs = []
        datasetRefs = butler.registry.queryDatasets("*_config", collections=collections)
        for datasetRef in datasetRefs:
            config = butler.getDirect(datasetRef)
            if isinstance(config, AnalysisBaseConfig):
                configs.append(config)
    plotNames = []
    for config in configs:
        connections: PipelineTaskConnections = config.connections.ConnectionsClass(config=config)
        for connection in iterConnections(connections, "outputs"):
            if connection.storageClass == "Plot":
                plotNames.append(connection.name)
    return plotNames
