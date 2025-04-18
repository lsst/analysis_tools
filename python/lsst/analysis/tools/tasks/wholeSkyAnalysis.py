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
    "WholeSkyAnalysisConfig",
    "WholeSkyAnalysisTask",
)
from typing import Mapping

from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask, AnalysisTool


class WholeSkyAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap",),
    defaultTemplates={"outputName": "objectTableCore_wholeSky", "inputName": "objectTableCore_metricsTable"},
):
    data = ct.Input(
        doc="Tract based object table to load from the butler",
        name="{inputName}",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("skymap",),
    )

    skymap = ct.Input(
        doc="The skymap that covers the tract that the data is from.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    def __init__(self, *, config: AnalysisBaseConfig = None):
        AnalysisBaseConnections.__init__(self, config=config)
        names: Mapping[str, AnalysisTool] = {}
        for action in config.atools:
            if action.dynamicOutputNames:
                outNames = action.getOutputNames(config=config)
            else:
                outNames = action.getOutputNames()
            names.update({name: action for name in outNames})

        for name, action in names.items():
            name = f"{self.config.connections.outputName}_{name}"
            # If the metric is band-specific, add band to dimensions.
            if any(f"_{band}_" in name for band in action.bands) and "band" not in self.dimensions:
                multiple = True
                dimensions = self.dimensions.union({"band"})
                for band in action.bands:
                    if f"_{band}_" in name:
                        name = name.replace(f"_{band}_", "_")
                        break
            else:
                multiple = False
                dimensions = self.dimensions

            outConnection = ct.Output(
                name=name,
                storageClass="Plot",
                doc="Dynamic connection for plotting",
                dimensions=dimensions,
                multiple=multiple,
            )
            setattr(self, name, outConnection)


class WholeSkyAnalysisConfig(AnalysisBaseConfig, pipelineConnections=WholeSkyAnalysisConnections):
    def setDefaults(self):
        super().setDefaults()
        self.bands = []


class WholeSkyAnalysisTask(AnalysisPipelineTask):
    ConfigClass = WholeSkyAnalysisConfig
    _DefaultName = "wholeSkyAnalysis"

    def run(self, *, data, **kwargs):
        data = dict(data)
        return super().run(data=data, **kwargs)
