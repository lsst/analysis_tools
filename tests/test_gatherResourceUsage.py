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

import os
from unittest import main

import lsst.utils.tests
from lsst.analysis.tools.tasks import GatherResourceUsageConfig, GatherResourceUsageTask
from lsst.daf.butler import Butler, DatasetType, DeferredDatasetHandle
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.pipe.base import TaskMetadata

TESTDIR = os.path.abspath(os.path.dirname(__file__))

test_metadata_json = """{"metadata": {"quantum":
            {"scalars":{"__version__":1},
             "arrays": {"endMaxResidentSetSize":[1234567890],
                        "prepCpuTime":[200.0],
                        "initCpuTime":[201.0],
                        "startCpuTime":[203.0],
                        "endCpuTime":[228.0],
                        "startUtc":["2025-01-08T22:20:00+00:00"],
                        "endUtc":["2025-01-08T22:20:45+00:00"]}}}}
                    """


class TestGatherResourceUsage(lsst.utils.tests.TestCase):
    """Tests for the GatherResourceUsage class and methods."""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)
        Butler.makeRepo(self.root)
        self.butler = Butler(self.root, run="test_run")

    def tearDown(self):
        removeTestTempDir(self.root)

    def test_gatherResourceUsage(self):
        test_dataset_type = DatasetType(
            "test_dataset_type", dimensions=[], universe=self.butler.dimensions, storageClass="TaskMetadata"
        )
        self.butler.registry.registerDatasetType(test_dataset_type)
        task_metadata = TaskMetadata.model_validate_json(test_metadata_json)
        ref = self.butler.put(task_metadata, test_dataset_type)
        config = GatherResourceUsageConfig(dimensions=[])
        gather_resource_usage_task = GatherResourceUsageTask(config=config)
        input_metadata = [
            DeferredDatasetHandle(butler=self.butler, ref=ref, storageClass="TaskMetadata", parameters=None)
        ]
        ru_output = gather_resource_usage_task.run(
            universe=self.butler.dimensions, input_metadata=input_metadata
        ).getDict()["output_table"]
        self.assertFloatsAlmostEqual(ru_output["memory"].values[0], 1.23456789e09, atol=1e-3, rtol=1e-4)
        self.assertFloatsAlmostEqual(ru_output["init_time"].values[0], 2.0, atol=1e-3, rtol=1e-4)
        self.assertFloatsAlmostEqual(ru_output["run_time"].values[0], 25.0, atol=1e-3, rtol=1e-4)
        self.assertFloatsAlmostEqual(ru_output["wall_time"].values[0], 45.0, atol=1e-3, rtol=1e-4)


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
