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
from collections.abc import Iterable
from typing import cast

from ..interfaces._interfaces import KeyedData
from lsst.daf.butler import DeferredDatasetHandle

from lsst.pipe.base.connections import InputQuantizedConnection, OutputQuantizedConnection

__all__ = (
    "DiaFakesVisitAnalysisConnections",
    "DiaFakesVisitAnalysisConfig",
    "DiaFakesVisitAnalysisTask",
)

from lsst.pipe.base import QuantumContext, connectionTypes as ct

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask

from astropy.table import Table, vstack


class DiaFakesVisitAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("visit", "band"),
    defaultTemplates={"coaddName": "goodSeeing", "fakesType": "fakes_"},
):
    data = ct.Input(
        doc="CcdVisit-based Matched fake to load from the butler",
        name="{fakesType}{coaddName}Diff_matchDiaSrc",
        storageClass="DataFrame",
        deferLoad=True,
        multiple=True,
        dimensions=("visit", "band", "detector"),
    )


class DiaFakesVisitAnalysisConfig(
    AnalysisBaseConfig, pipelineConnections=DiaFakesVisitAnalysisConnections
):
    pass


class DiaFakesVisitAnalysisTask(AnalysisPipelineTask):
    ConfigClass = DiaFakesVisitAnalysisConfig
    _DefaultName = "FakesDiaVisitAnalysis"

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection
    ) -> None:

        inputs = butlerQC.get(inputRefs)
        fakesTable = self.loadMultipleData(inputs["data"])
        inputs["data"] = fakesTable
        outputs = self.run(**inputs)

        butlerQC.put(outputs, outputRefs)

    def loadMultipleData(
        self,
        inputs: Iterable[DeferredDatasetHandle],
        names: Iterable[str] | None = None
    ) -> KeyedData:

        fakesTables = []
        for aninput in inputs:
            fakesTables.append(Table.from_pandas(aninput.get()))

        fakesTable = vstack(fakesTables)

        return cast(KeyedData, fakesTable.to_pandas())
