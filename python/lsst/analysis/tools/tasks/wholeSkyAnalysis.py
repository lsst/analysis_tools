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

from lsst.pex.config import ListField
from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class WholeSkyAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap",),
    defaultTemplates={"outputName": "objectTableCore_wholeSky", "inputName": "objectTableCore_metricsTable"},
):
    data = ct.Input(
        doc="Tract based table to load from the butler",
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

    def __init__(self, *, config=None):
        """Customize the dimensions of the inputs/outputs for a specific
        instance. This enables it to be dynamically set at runtime,
        allowing the task to work with different datasets.

        Parameters
        ----------
        config : `WholeSkyAnalysisConfig`
            A config for `WholeSkyAnalysisConfig`.
        """
        super().__init__(config=config)
        if config.inputDimensions:
            self.data = ct.Input(
                name=self.data.name,
                doc=self.data.doc,
                storageClass=self.data.storageClass,
                dimensions=frozenset(sorted(config.outputDimensions)),
                deferLoad=self.data.deferLoad,
                multiple=self.data.multiple,
            )

        if config.outputDimensions:
            self.dimensions.clear()
            self.dimensions.update(frozenset(sorted(config.outputDimensions)))


class WholeSkyAnalysisConfig(AnalysisBaseConfig, pipelineConnections=WholeSkyAnalysisConnections):
    inputDimensions = ListField[str](
        default=(),
        doc=("Override the dimensions of the input data. "),
    )
    outputDimensions = ListField[str](
        default=(), doc=("Override the dimensions of the output data." "Also overrides the task dimensions.")
    )


class WholeSkyAnalysisTask(AnalysisPipelineTask):
    ConfigClass = WholeSkyAnalysisConfig
    _DefaultName = "wholeSkyAnalysis"
