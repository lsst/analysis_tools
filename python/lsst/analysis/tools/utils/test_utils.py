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


from lsst.ip.isr.isrTask import IsrTask
from lsst.pipe.base import Pipeline
from lsst.pipe.base.pipelineIR import LabeledSubset


def make_test_reference_pipeline():
    """Make a test reference pipeline containing initial single-frame tasks."""
    reference_pipeline = Pipeline("reference_pipeline")
    reference_pipeline.addTask(IsrTask, "isr")
    reference_pipeline._pipelineIR.labeled_subsets["test_subset"] = LabeledSubset("test_subset", set(), None)
    reference_pipeline.addLabelToSubset("test_subset", "isr")
    return reference_pipeline
