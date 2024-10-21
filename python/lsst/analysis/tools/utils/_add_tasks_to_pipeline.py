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

__all__ = ["add_tasks_to_pipeline"]

import copy
import logging

from lsst.pipe.base import Pipeline


def add_tasks_to_pipeline(
    reference_pipeline: Pipeline | str,
    input_pipelines: list[Pipeline] | list[str],
    subset_name: str | None = None,
    new_subset_description: str = "",
    instrument: str | None = None,
    log_level: int = logging.INFO,
) -> Pipeline:
    """Add task(s) to a reference pipeline.

    Parameters
    ----------
    reference_pipeline: Pipeline | `str`
        Location of a reference pipeline definition YAML file.
    input_pipelines: `list`[Pipeline] | `list`[`str`]
        Location(s) of input pipeline definition YAML file(s). Tasks from
        input_pipelines will be added to reference_pipeline.
    subset_name: `str`, optional
        All tasks from input_pipelines will be added to this subset. If the
        subset does not exist it will be created.
    new_subset_description: `str`, optional
        The description for the new subset.
    instrument : `str`, optional
        Add instrument overrides. Must be a fully qualified class name.
    log_level : `int`, optional
        The log level to use for logging.

    Returns
    -------
    pipeline : `lsst.pipe.base.Pipeline`
        An expanded pipeline with the tasks and subsets from reference_pipeline
        and input_pipelines.
    """
    # Instantiate logger.
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    # Load the pipeline and apply config overrides, if supplied.
    # Copy to preserve the original pipeline objects.
    if isinstance(reference_pipeline, str):
        pipeline = Pipeline.from_uri(reference_pipeline)
    else:
        pipeline = copy.deepcopy(reference_pipeline)

    # Add an instrument override, if provided.
    if instrument:
        pipeline.addInstrument(instrument)

    # Record all input task labels and merge all input pipelines into the
    # reference pipeline.
    input_tasks = set()
    for input_pipeline in input_pipelines:
        if isinstance(input_pipeline, str):
            input_pipeline = Pipeline.fromFile(input_pipeline)
        input_tasks.update(input_pipeline.task_labels)
        pipeline.mergePipeline(input_pipeline)

    # Add all tasks to subset_name. If the subset does not exist, create it.
    if isinstance(subset_name, str):
        if subset_name in pipeline.subsets.keys():
            for input_task in input_tasks:
                pipeline.addLabelToSubset(subset_name, input_task)
                subset_grammar = f"the existing subset {subset_name}"
        else:
            pipeline.addLabeledSubset(subset_name, new_subset_description, input_tasks)
            subset_grammar = f"a new subset {subset_name}"

        # Logging info.
        task_grammar = "task" if len(input_tasks) == 1 else "tasks"
        logger.info(
            "Added %d %s to %s",
            len(input_tasks),
            task_grammar,
            subset_grammar,
        )

    return pipeline
