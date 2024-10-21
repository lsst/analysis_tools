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

import logging
import os
import unittest

import lsst.utils.tests
from lsst.analysis.tools import add_tasks_to_pipeline
from lsst.analysis.tools.utils.test_utils import make_test_reference_pipeline
from lsst.pipe.base import Pipeline
from lsst.pipe.tasks.calibrate import CalibrateTask
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
from lsst.utils.tests import MemoryTestCase, TestCase

TEST_DIR = os.path.abspath(os.path.dirname(__file__))


class PipelineToolsUtilsTestCase(TestCase):
    """Test the utility functions in the source_injection package."""

    def setUp(self):
        self.reference_pipeline = make_test_reference_pipeline()

    def tearDown(self):
        del self.reference_pipeline

    def test_add_tasks_to_pipeline(self):
        input_pipeline1 = Pipeline("input_pipeline1")
        input_pipeline1.addTask(CharacterizeImageTask, "characterizeImage")
        input_pipeline2 = Pipeline("input_pipeline2")
        input_pipeline2.addTask(CalibrateTask, "calibrate")

        # Merge the injection pipeline into the main reference pipeline.
        new_subset_name = "new_test_subset"
        new_pipeline = add_tasks_to_pipeline(
            reference_pipeline=self.reference_pipeline,
            input_pipelines=[input_pipeline1, input_pipeline2],
            subset_name=new_subset_name,
            new_subset_description="lorem ipsum",
            instrument="lsst.obs.subaru.HyperSuprimeCam",
            log_level=logging.DEBUG,
        )

        # Test that all expected labels are present in the new pipeline.
        reference_task_labels = self.reference_pipeline.task_labels
        expected_task_labels = (
            set(reference_task_labels) | set(input_pipeline1.task_labels) | set(input_pipeline2.task_labels)
        )
        actual_task_labels = set(new_pipeline.task_labels)
        self.assertEqual(expected_task_labels, actual_task_labels)

        # Test that the new pipeline has a new subset with the correct tasks.
        expected_new_subset_labels = set(input_pipeline1.task_labels) | set(input_pipeline2.task_labels)
        actual_new_subset_labels = new_pipeline.subsets[new_subset_name]
        self.assertEqual(expected_new_subset_labels, actual_new_subset_labels)

        # Test that tasks can be added to an existing subset.
        # Grab the name of the first subset in the reference pipeline.
        subset_name = next(iter(self.reference_pipeline.subsets))
        new_pipeline = add_tasks_to_pipeline(
            reference_pipeline=self.reference_pipeline,
            input_pipelines=[input_pipeline1, input_pipeline2],
            subset_name=subset_name,
            new_subset_description="",
            instrument="lsst.obs.subaru.HyperSuprimeCam",
            log_level=logging.DEBUG,
        )
        expected_subset_labels = (
            self.reference_pipeline.subsets[subset_name]
            | set(input_pipeline1.task_labels)
            | set(input_pipeline2.task_labels)
        )
        actual_subset_labels = new_pipeline.subsets[subset_name]
        self.assertEqual(expected_subset_labels, actual_subset_labels)


class MemoryTestCase(MemoryTestCase):
    """Test memory usage of functions in this script."""

    pass


def setup_module(module):
    """Configure pytest."""
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
