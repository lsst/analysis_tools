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

__all__ = ("ObjectTableSurveyAnalysisTask",)


from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from lsst.daf.butler import DataCoordinate, DeferredDatasetHandle

from astropy.table import vstack

from lsst.pipe.base import connectionTypes as ct
from lsst.skymap import BaseSkyMap

from ..actions.plot.plotUtils import shorten_list
from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask, KeyedData


class ObjectTableSurveyAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("skymap",),
    defaultTemplates={"input": "deepCoadd"},
):
    skymap = ct.Input(
        doc="Input definition of geometry/bbox and projection/wcs for warped exposures",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    data = ct.Input(
        doc="Input catalog of objects",
        name="objectTable_tract",
        storageClass="ArrowAstropy",
        deferLoad=True,
        dimensions=("tract", "skymap"),
        multiple=True,
    )


class ObjectTableSurveyAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=ObjectTableSurveyAnalysisConnections
):
    pass


class ObjectTableSurveyAnalysisTask(AnalysisPipelineTask):
    """A specialized ``AnalysisPipelineTask`` for multiple tracts."""

    ConfigClass = ObjectTableSurveyAnalysisConfig
    _DefaultName = "objectTableSurveyAnalysis"

    def parsePlotInfo(
        self, inputs: Mapping[str, Any] | None, dataId: DataCoordinate | None, connectionName: str = "data"
    ) -> Mapping[str, str]:
        # Docstring inherited

        if inputs is None:
            tableName = ""
            run = ""
            tracts = []
        else:
            tableName = inputs[connectionName][0].ref.datasetType.name
            run = inputs[connectionName][0].ref.run
            tracts = [data.ref.dataId["tract"] for data in list(inputs[connectionName])]

        # Initialize the plot info dictionary
        plotInfo = {
            "tableName": tableName,
            "run": run,
            "tract": shorten_list(tracts),
        }

        self._populatePlotInfoWithDataId(plotInfo, dataId)
        return plotInfo

    def loadData(  # type: ignore[override]
        self,
        handle: Iterable[DeferredDatasetHandle],
        names: Iterable[str] | None = None,
    ) -> KeyedData:
        """Load the minimal set of keyed data from the input dataset.

        Parameters
        ----------
        handle : `Iterable` of `DeferredDatasetHandle`
            Handle to load the dataset with only the specified columns.
        names : `Iterable` of `str`
            The names of keys to extract from the dataset.
            If `names` is `None` then the `collectInputNames` method
            is called to generate the names.
            For most purposes these are the names of columns to load from
            a catalog or data frame.

        Returns
        -------
        result: `KeyedData`
            The dataset with only the specified keys loaded.
        """
        # Except associatedSourceTractAnalysis, all other tasks trivially
        # subclass this, so names don't get utilized.
        if names is None:
            names = self.collectInputNames()

        cats = []
        for h in handle:
            cats.append(h.get(parameters={"columns": names}))
        return vstack(cats)
