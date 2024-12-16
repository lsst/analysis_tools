# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__all__ = (
    "ConsolidateResourceUsageConfig",
    "ConsolidateResourceUsageConnections",
    "ConsolidateResourceUsageTask",
    "GatherResourceUsageConfig",
    "GatherResourceUsageConnections",
    "GatherResourceUsageTask",
    "ResourceUsageQuantumGraphBuilder",
)

import argparse
import dataclasses
import datetime
import logging
import re
from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import pandas as pd
from astropy.time import Time
from lsst.daf.butler import Butler, DatasetRef, DatasetType
from lsst.pex.config import Field, ListField
from lsst.pipe.base import (
    Instrument,
    PipelineTask,
    PipelineTaskConfig,
    PipelineTaskConnections,
    QuantumGraph,
    Struct,
)
from lsst.pipe.base import connectionTypes as cT
from lsst.pipe.base.pipeline_graph import PipelineGraph
from lsst.pipe.base.quantum_graph_builder import QuantumGraphBuilder
from lsst.pipe.base.quantum_graph_skeleton import DatasetKey, QuantumGraphSkeleton

# It's not great to be importing a private symbol, but this is a temporary
# workaround for the fact that prior to w.2022.10, the units for memory values
# written in task metadata were platform-dependent.  Once we no longer care
# about older runs, this import and the code that uses it can be removed.
from lsst.utils.usage import _RUSAGE_MEMORY_MULTIPLIER

_LOG = logging.getLogger(__name__)


class ConsolidateResourceUsageConnections(PipelineTaskConnections, dimensions=()):
    """Connection definitions for `ConsolidateResourceUsageTask`."""

    output_table = cT.Output(
        name="ResourceUsageSummary",
        storageClass="DataFrame",
        dimensions=(),
        doc="Consolidated table of resource usage statistics. One row per task label",
    )

    def __init__(self, *, config):
        super().__init__(config=config)
        for name in self.config.input_names:
            setattr(
                self,
                name,
                cT.Input(
                    name,
                    storageClass="DataFrame",
                    dimensions=(),
                    doc="Resource usage statistics for a task.",
                ),
            )
            self.inputs.add(name)


class ConsolidateResourceUsageConfig(
    PipelineTaskConfig, pipelineConnections=ConsolidateResourceUsageConnections
):
    """Configuration definitions for `ConsolidateResourceUsageTask`."""

    input_names = ListField[str](
        doc="Input resource usage dataset type names",
        default=[],
    )


class ConsolidateResourceUsageTask(PipelineTask):
    """A `PipelineTask` that summarizes task resource usage into a single
    table with per-task rows.

    Notes
    -----
    This is an unusual `PipelineTask` in that its input connection has
    dynamic dimensions, and its quanta are generally built via a custom
    quantum-graph builder defined in the same module.
    """

    ConfigClass = ConsolidateResourceUsageConfig
    _DefaultName = "consolidateResourceUsage"

    def run(self, **kwargs: Any) -> Struct:
        quantiles = []
        for input_name, ru_table in kwargs.items():
            if not input_name.endswith("resource_usage"):
                continue
            else:
                df = ru_table.quantile(
                    [0.0, 0.01, 0.05, 0.32, 0.50, 0.68, 0.95, 0.99, 1.0],
                    numeric_only=True,
                ).reset_index()
                df["task"] = input_name.replace("_resource_usage", "")
                df["quanta"] = len(ru_table)
                df["integrated_runtime"] = ru_table["run_time"].sum()

                quantiles.append(
                    df[
                        [
                            "index",
                            "quanta",
                            "task",
                            "memory",
                            "init_time",
                            "run_time",
                            "wall_time",
                            "integrated_runtime",
                        ]
                    ]
                )

        full_quantiles = pd.concat(quantiles)
        full_quantiles["percentile"] = (full_quantiles["index"] * 100).astype(int)
        full_quantiles["percentile_name"] = "p" + full_quantiles["percentile"].astype(str).str.zfill(3)
        full_quantiles["memoryGB"] = full_quantiles["memory"] / 1024 / 1024 / 1024
        full_quantiles["integrated_runtime_hrs"] = full_quantiles["integrated_runtime"] / 3600.0
        memoryGB = pd.pivot_table(
            full_quantiles, values="memoryGB", columns=["percentile_name"], index=["task"]
        ).add_prefix("mem_GB_")
        runtime = pd.pivot_table(
            full_quantiles, values="run_time", columns=["percentile_name"], index=["task"]
        ).add_prefix("runtime_s_")
        walltime = pd.pivot_table(
            full_quantiles, values="wall_time", columns=["percentile_name"], index=["task"]
        ).add_prefix("walltime_s_")
        memrun = pd.merge(
            pd.merge(
                memoryGB.reset_index(),
                runtime.reset_index(),
                left_on="task",
                right_on="task",
            ),
            walltime.reset_index(),
            left_on="task",
            right_on="task",
        )
        memrun = pd.merge(
            full_quantiles[["task", "quanta", "integrated_runtime_hrs"]]
            .drop_duplicates()
            .sort_values("task"),
            memrun,
        )

        return Struct(output_table=memrun)


class GatherResourceUsageConnections(
    PipelineTaskConnections, dimensions=(), defaultTemplates={"input_task_label": "PLACEHOLDER"}
):
    """Connection definitions for `GatherResourceUsageTask`."""

    output_table = cT.Output(
        "{input_task_label}_resource_statistics",  # Should always be overridden.
        storageClass="DataFrame",
        dimensions=(),
        doc=(
            "Table that aggregates memory and CPU usage statistics from one "
            "or more tasks. "
            "This will have one row for each data ID, with columns for each "
            "task or method's memory usage and runtime."
        ),
    )
    input_metadata = cT.Input(
        "{input_task_label}_metadata",  # Should always be overridden.
        storageClass="TaskMetadata",
        dimensions=(),  # Actually set in __init__, according to configuration.
        doc="Metadata dataset for another task to gather resource usage from.",
        multiple=True,
        deferLoad=True,
    )

    def __init__(self, *, config):
        super().__init__(config=config)
        if "PLACEHOLDER" in self.output_table.name:
            raise ValueError("Connection configuration for output_table must be overridden.")
        if "PLACEHOLDER" in self.input_metadata.name:
            raise ValueError("Connection configuration for input_metadata must be overridden.")
        # Override the empty dimension set the connection was defined with with
        # those the task was configured with.
        self.input_metadata = dataclasses.replace(
            self.input_metadata,
            dimensions=list(self.config.dimensions),
        )


class GatherResourceUsageConfig(PipelineTaskConfig, pipelineConnections=GatherResourceUsageConnections):
    """Configuration definitions for `GatherResourceUsageTask`."""

    dimensions = ListField[str](
        doc=(
            "The quantum dimensions for the input metadata connection, and "
            "the columns (after expansion to include implied dimensions) used "
            "to identify rows in the output table."
        ),
    )
    memory = Field[bool](
        doc=(
            "Whether to extract peak memory usage (maximum resident set size) "
            "for this task. "
            "Note that memory usage cannot be further subdivided because only "
            "a per-process peak is available (and hence if multiple quanta "
            "are run in one quantum, even per-quantum values may be "
            "misleading)."
        ),
        default=True,
    )
    prep_time = Field[bool](
        doc=(
            "Whether to extract the CPU time duration for the work the "
            "middleware does prior to initializing the task (mostly checking "
            "for input dataset existence)."
        ),
        default=False,
    )
    init_time = Field[bool](
        doc=("Whether to extract the CPU time duration for actually constructing the task."),
        default=True,
    )
    run_time = Field[bool](
        doc=("Whether to extract the CPU time duration for actually executing the task."),
        default=True,
    )
    wall_time = Field[bool](
        doc=("Whether to extract the wall_time duration for actually executing the task."),
        default=True,
    )
    method_times = ListField[str](
        doc=(
            "Names of @lsst.utils.timer.timeMethod-decorated methods for "
            "which CPU time durations should also be extracted.  Use '.' "
            "separators to refer to subtask methods at arbitrary depth."
        ),
        optional=False,
        default=[],
    )
    input_task_label = Field[str](
        doc=(
            "Label for the top-level task whose metadata is being processed "
            "within its own metadata file, if this differs from the prefix of "
            "connections.input_metadata."
        ),
        default=None,
        optional=True,
    )


class GatherResourceUsageTask(PipelineTask):
    """A `PipelineTask` that gathers resource usage statistics from task
    metadata.

    Notes
    -----
    This is an unusual `PipelineTask` in that its input connection has
    dynamic dimensions.

    Its output table has columns for each of the dimensions of the input
    metadata's data ID, as well as (subject to configuration):

    - ``memory``: the maximum resident set size for the entire quantum
      (in bytes);
    - ``prep_time``: the time spent in the pre-initialization step in
      which the middleware checks which of the quantum's inputs are available;
    - ``init_time``: the time spent in task construction;
    - ``run_time``: the time spent executing the task's runQuantum
      method.
    - ``wall_time`` : elapsed time in the pre-initialization step, in task
      construction, and in executing the task's runQuantum method.
      Specifically, this is the difference between `prepUtc`, which triggers
      as soon as single quantum execution has begun (but can include some
      checks and running `updatedQuantumInputs`), and `endUtc`, which triggers
      immediately after `runQuantum`.
    - ``{method}``: the time spent in a particular task or subtask
      method decorated with `lsst.utils.timer.timeMethod`.

    All time durations are CPU times in seconds, and all columns are 64-bit
    floating point.  Methods or steps that did not run are given a duration of
    zero.

    It is expected that this task will be configured to run multiple times in
    most pipelines, often once for each other task in the pipeline.
    """

    ConfigClass = GatherResourceUsageConfig
    _DefaultName = "gatherResourceUsage"

    def runQuantum(
        self,
        butlerQC,
        inputRefs,
        outputRefs,
    ):
        # Docstring inherited.
        # This override exists just so we can pass the butler registry's
        # DimensionUniverse to run in order to standardize the dimensions.
        inputs = butlerQC.get(inputRefs)
        outputs = self.run(butlerQC.dimensions, **inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, universe, input_metadata):
        """Gather resource usage statistics from per-quantum metadata.

        Parameters
        ----------
        universe : `DimensionUniverse`
            Object managing all dimensions recognized by the butler; used to
            standardize and expand `GatherResourceUsageConfig.dimensions`.
        input_metadata : `list` [ `DeferredDatasetHandle` ]
            List of `lsst.daf.butler.DeferredDatasetHandle` that can be used to
            load all input metadata datasets.

        Returns
        -------
        result : `Struct`
            Structure with a single element:

            - ``outout_table``: a `pandas.DataFrame` that aggregates the
              configured resource usage statistics.
        """
        dimensions = universe.conform(self.config.dimensions)
        # Transform input list into a dict keyed by data ID.
        handles_by_data_id = {}
        for handle in input_metadata:
            handles_by_data_id[handle.dataId] = handle
        n_rows = len(handles_by_data_id)
        # Create a dict of empty column arrays that we'll ultimately make into
        # a table.
        columns = {
            d: np.zeros(n_rows, dtype=_dtype_from_field_spec(universe.dimensions[d].primaryKey))
            for d in dimensions.names
        }
        for attr_name in ("memory", "prep_time", "init_time", "run_time", "wall_time"):
            if getattr(self.config, attr_name):
                columns[attr_name] = np.zeros(n_rows, dtype=float)
        for method_name in self.config.method_times:
            columns[method_name] = np.zeros(n_rows, dtype=float)
        # Populate the table, one row at a time.
        warned_about_metadata_version = False
        for index, (data_id, handle) in enumerate(handles_by_data_id.items()):
            # Fill in the data ID columns.
            for k, v in data_id.mapping.items():
                columns[k][index] = v
            # Load the metadata dataset and fill in the columns derived from
            # it.
            metadata = handle.get()
            try:
                quantum_metadata = metadata["quantum"]
            except KeyError:
                self.log.warning(
                    "Metadata dataset %s @ %s has no 'quantum' key.",
                    handle.ref.datasetType.name,
                    handle.dataId,
                )
            else:
                if self.config.memory:
                    columns["memory"][index], warned_about_metadata_version = self._extract_memory(
                        quantum_metadata,
                        handle,
                        warned_about_metadata_version,
                    )
                for key, value in self._extract_quantum_timing(quantum_metadata).items():
                    columns[key][index] = value
            for key, value in self._extract_method_timing(metadata, handle).items():
                columns[key][index] = value
        return Struct(output_table=pd.DataFrame(columns, copy=False))

    def _extract_memory(self, quantum_metadata, handle, warned_about_metadata_version):
        """Extract maximum memory usage from quantum metadata.

        Parameters
        ----------
        quantum_metadata : `lsst.pipe.base.TaskMetadata`
            The nested metadata associated with the label "quantum" inside a
            PipelineTask's metadata.
        handle : `lsst.daf.butler.DeferredDatasetHandle`
            Butler handle for the metadata dataset; used to identify the
            metadata in diagnostic messages only.
        warned_about_metadata_version : `bool`
            Whether we have already emitted at least one warning about old
            metadata versions.

        Returns
        -------
        memory : `float`
            Maximum memory usage in bytes.
        warned_about_metadata_version : `bool`
            Whether we have now emitted at least one warning about old
            metadata versions.
        """
        # Attempt to work around memory units being
        # platform-dependent for metadata written prior to
        # w.2022.10.
        memory_multiplier = 1
        if quantum_metadata.get("__version__", 0) < 1:
            memory_multiplier = _RUSAGE_MEMORY_MULTIPLIER
            msg = (
                "Metadata dataset %s @ %s is too old; guessing memory units by "
                "assuming the platform has not changed"
            )
            if not warned_about_metadata_version:
                self.log.warning(msg, handle.ref.datasetType.name, handle.dataId)
                self.log.warning(
                    "Warnings about memory units for other inputs " "will be emitted only at DEBUG level."
                )
                warned_about_metadata_version = True
            else:
                self.log.debug(msg, handle.ref.datasetType.name, handle.dataId)
        return (
            quantum_metadata["endMaxResidentSetSize"] * memory_multiplier,
            warned_about_metadata_version,
        )

    def _extract_quantum_timing(self, quantum_metadata):
        """Extract timing for standard PipelineTask quantum-execution steps
        from metadata.

        Parameters
        ----------
        quantum_metadata : `lsst.pipe.base.TaskMetadata`
            The nested metadata associated with the label "quantum" inside a
            PipelineTask's metadata.

        Returns
        -------
        timing : `dict` [ `str`, `float` ]
            CPU times in bytes, for all stages enabled in configuration.
        """
        end_time = quantum_metadata["endCpuTime"]
        times = [
            quantum_metadata["prepCpuTime"],
            quantum_metadata.get("initCpuTime", end_time),
            quantum_metadata.get("startCpuTime", end_time),
            end_time,
        ]

        quantum_timing = {
            attr_name: end - begin
            for attr_name, begin, end in zip(
                ["prep_time", "init_time", "run_time"],
                times[:-1],
                times[1:],
            )
            if getattr(self.config, attr_name)
        }
        if self.config.wall_time:
            start_wall_time = Time(quantum_metadata["prepUtc"].split("+")[0])
            end_wall_time = Time(quantum_metadata["endUtc"].split("+")[0])
            quantum_timing["wall_time"] = (end_wall_time - start_wall_time).sec

        return quantum_timing

    def _extract_method_timing(self, metadata, handle):
        """Extract timing for standard PipelineTask quantum-execution steps
        from metadata.

        Parameters
        ----------
        quantum_metadata : `lsst.pipe.base.TaskMetadata`
            The nested metadata associated with the label "quantum" inside a
            PipelineTask's metadata.
        handle : `lsst.daf.butler.DeferredDatasetHandle`
            Butler handle for the metadata dataset; used infer the prefix used
            for method names within the metadata.

        Returns
        -------
        timing : `dict` [ `str`, `float` ]
            CPU times in bytes, for all methods enabled in configuration.
        """
        if self.config.input_task_label is not None:
            task_label = self.config.input_task_label
        else:
            task_label = handle.ref.datasetType.name[: -len("_metadata")]
        result = {}
        for method_name in self.config.method_times:
            terms = [task_label] + list(method_name.split("."))
            metadata_method_name = ":".join(terms[:-1]) + "." + terms[-1]
            try:
                method_start_time = metadata[f"{metadata_method_name}StartCpuTime"]
                method_end_time = metadata[f"{metadata_method_name}EndCpuTime"]
            except KeyError:
                # A method missing from the metadata is not a problem;
                # it's reasonable for configuration or even runtime
                # logic to result in a method not being called.  When
                # that happens, we just let the times stay zero.
                pass
            else:
                result[f"{task_label}.{method_name}"] = method_end_time - method_start_time
        return result


def _dtype_from_field_spec(field_spec):
    """Return the `np.dtype` that can be used to hold the values of a butler
    dimension field.

    Parameters
    ----------
    field_spec : `lsst.daf.butler.core.ddl.FieldSpec`
        Object describing the field in a SQL-friendly sense.

    Returns
    -------
    dtype : `np.dtype`
        Numpy data type description.
    """
    python_type = field_spec.getPythonType()
    if python_type is str:
        return np.dtype((str, field_spec.length))
    else:
        return np.dtype(python_type)


class ResourceUsageQuantumGraphBuilder(QuantumGraphBuilder):
    """Custom quantum graph generator and pipeline builder for resource
    usage summary tasks.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler client to query for inputs and dataset types.
    dataset_type_names : `~collections.abc.Iterable` [ `str` ], optional
        Iterable of dataset type names or shell-style glob patterns for the
        metadata datasets to be used as input.  Default is all datasets ending
        with ``_metadata`` (other than the resource-usage summary tasks' own
        metadata outputs, where are always ignored).  A gather-resource task
        with a single quantum is created for each matching metadata dataset.
    where : `str`, optional
        Data ID expression that constrains the input metadata datasets.
    input_collections : `~collections.abc.Sequence` [ `str` ], optional
        Sequence of collections to search for inputs.  If not provided,
        ``butler.collections`` is used and must not be empty.
    output_run : `str`, optional
        Output `~lsst.daf.butler.CollectionType.RUN` collection name.  If not
        provided, ``butler.run`` is used and must not be `None`.
    skip_existing_in : `~collections.abc.Sequence` [ `str` ], optional
        Sequence of collections to search for outputs, allowing quanta whose
        outputs exist to be skipped.
    clobber : `bool`, optional
        Whether *execution* of this quantum graph will permit clobbering.  If
        `False` (default), existing outputs in ``output_run`` are an error
        unless ``skip_existing_in`` will cause those quanta to be skipped.

    Notes
    -----
    The resource usage summary tasks cannot easily be added to a regular
    pipeline, as it's much more natural to have the gather tasks run
    automatically on all *other* tasks. And we can generate a quantum graph
    for these particular tasks much more efficiently than the general-purpose
    algorithm could.
    """

    def __init__(
        self,
        butler: Butler,
        *,
        dataset_type_names: Iterable[str] | None = None,
        where: str = "",
        input_collections: Sequence[str] | None = None,
        output_run: str | None = None,
        skip_existing_in: Sequence[str] = (),
        clobber: bool = False,
    ):
        # Start by querying for metadata datasets, since we'll need to know
        # which dataset types exist in the input collections in order to
        # build the pipeline.
        input_dataset_types: Any
        if not dataset_type_names:
            input_dataset_types = "*_metadata"
        else:
            input_dataset_types = dataset_type_names
        pipeline_graph = PipelineGraph()
        metadata_refs: dict[str, set[DatasetRef]] = {}
        consolidate_config = ConsolidateResourceUsageConfig()
        for results in butler.registry.queryDatasets(
            input_dataset_types,
            where=where,
            findFirst=True,
            collections=input_collections,
        ).byParentDatasetType():
            input_metadata_dataset_type = results.parentDatasetType
            refs_for_type = set(results)
            if refs_for_type:
                gather_task_label, gather_dataset_type_name = self._add_gather_task(
                    pipeline_graph, input_metadata_dataset_type
                )
                metadata_refs[gather_task_label] = refs_for_type
                consolidate_config.input_names.append(gather_dataset_type_name)
        pipeline_graph.add_task(
            task_class=ConsolidateResourceUsageTask,
            config=consolidate_config,
            label=ConsolidateResourceUsageTask._DefaultName,
        )
        # Now that we have the pipeline graph, we can delegate to super.
        super().__init__(
            pipeline_graph,
            butler,
            input_collections=input_collections,
            output_run=output_run,
            skip_existing_in=skip_existing_in,
            clobber=clobber,
        )
        # We've already queried for all of our input datasets, so we don't want
        # to do that again in process_subgraph, even though that's where most
        # QG builders do their queries.
        self.gather_inputs: dict[str, list[DatasetKey]] = {}
        for gather_task_label, gather_input_refs in metadata_refs.items():
            gather_inputs_for_task: list[DatasetKey] = []
            for ref in gather_input_refs:
                dataset_key = DatasetKey(ref.datasetType.name, ref.dataId.required_values)
                self.existing_datasets.inputs[dataset_key] = ref
                gather_inputs_for_task.append(dataset_key)
            self.gather_inputs[gather_task_label] = gather_inputs_for_task

    @classmethod
    def _add_gather_task(
        cls, pipeline_graph: PipelineGraph, input_metadata_dataset_type: DatasetType
    ) -> tuple[str, str]:
        """Add a single configuration of `GatherResourceUsageTask` to a
        pipeline graph.

        Parameters
        ----------
        pipeline_graph : `lsst.pipe.base.PipelineGraph`
            Pipeline graph to modify in-place.
        input_metadata_dataset_type : `lsst.daf.butler.DatasetType`
            Dataset type for the task's input dataset, which is the metadata
            output of the task whose resource usage information is being
            extracted.

        Returns
        -------
        gather_task_label : `str`
            Label of the new task in the pipeline.
        gather_dataset_type_name : `str
            Name of the task's output table dataset type.
        """
        if (m := re.fullmatch(r"^(\w+)_metadata$", input_metadata_dataset_type.name)) is None:
            return
        elif "gatherResourceUsage" in input_metadata_dataset_type.name:
            return
        else:
            input_task_label = m.group(1)
        gather_task_label = f"{input_task_label}_gatherResourceUsage"
        gather_dataset_type_name = f"{input_task_label}_resource_usage"
        gather_config = GatherResourceUsageConfig()
        gather_config.dimensions = input_metadata_dataset_type.dimensions.names
        gather_config.connections.input_metadata = input_metadata_dataset_type.name
        gather_config.connections.output_table = gather_dataset_type_name
        pipeline_graph.add_task(
            label=gather_task_label,
            task_class=GatherResourceUsageTask,
            config=gather_config,
        )
        return gather_task_label, gather_dataset_type_name

    def process_subgraph(self, subgraph: PipelineGraph) -> QuantumGraphSkeleton:
        skeleton = QuantumGraphSkeleton(subgraph.tasks.keys())
        consolidate_inputs = []
        for task_node in subgraph.tasks.values():
            if task_node.task_class is GatherResourceUsageTask:
                quantum_key = skeleton.add_quantum_node(task_node.label, self.empty_data_id)
                skeleton.add_input_edges(quantum_key, self.gather_inputs[task_node.label])
                for write_edge in task_node.iter_all_outputs():
                    output_node = subgraph.dataset_types[write_edge.parent_dataset_type_name]
                    assert (
                        output_node.dimensions == self.universe.empty
                    ), "All outputs should have empty dimensions."
                    gather_output_key = skeleton.add_dataset_node(
                        write_edge.parent_dataset_type_name, self.empty_data_id
                    )
                    skeleton.add_output_edge(quantum_key, gather_output_key)
                    if write_edge.connection_name in task_node.outputs:
                        # Not a special output like metadata or log.
                        consolidate_inputs.append(gather_output_key)
            else:
                assert task_node.task_class is ConsolidateResourceUsageTask
                quantum_key = skeleton.add_quantum_node(task_node.label, self.empty_data_id)
                skeleton.add_input_edges(quantum_key, consolidate_inputs)
                for write_edge in task_node.iter_all_outputs():
                    output_node = subgraph.dataset_types[write_edge.parent_dataset_type_name]
                    assert (
                        output_node.dimensions == self.universe.empty
                    ), "All outputs should have empty dimensions."
                    consolidate_output_key = skeleton.add_dataset_node(
                        write_edge.parent_dataset_type_name, self.empty_data_id
                    )
                    skeleton.add_output_edge(quantum_key, consolidate_output_key)
        # We don't need to do any follow-up searches for output datasets,
        # because the outputs all have empty dimensions and the base
        # QuantumGraphBuilder takes care of those.
        return skeleton

    @classmethod
    def make_argument_parser(cls) -> argparse.ArgumentParser:
        """Make the argument parser for the command-line interface."""
        parser = argparse.ArgumentParser(
            description=(
                "Build a QuantumGraph that gathers and consolidates "
                "resource usage tables from existing metadata datasets."
            ),
        )
        parser.add_argument("repo", type=str, help="Path to data repository or butler configuration.")
        parser.add_argument("filename", type=str, help="Output filename for QuantumGraph.")
        parser.add_argument(
            "collections",
            type=str,
            nargs="+",
            help="Collection(s)s to search for input metadata.",
        )
        parser.add_argument(
            "--dataset-types",
            type=str,
            action="extend",
            help="Glob-style patterns for input metadata dataset types.",
        )
        parser.add_argument(
            "--where",
            type=str,
            default="",
            help="Data ID expression used when querying for input metadata datasets.",
        )
        parser.add_argument(
            "--output",
            type=str,
            help=(
                "Name of the output CHAINED collection. If this options is specified and "
                "--output-run is not, then a new RUN collection will be created by appending "
                "a timestamp to the value of this option."
            ),
            default=None,
            metavar="COLL",
        )
        parser.add_argument(
            "--output-run",
            type=str,
            help=(
                "Output RUN collection to write resulting images. If not provided "
                "then --output must be provided and a new RUN collection will be created "
                "by appending a timestamp to the value passed with --output."
            ),
            default=None,
            metavar="RUN",
        )
        return parser

    @classmethod
    def main(cls) -> None:
        """Run the command-line interface for this quantum-graph builder.

        This function provides the implementation for the
        ``build-gather-resource-usage-qg`` script.
        """
        parser = cls.make_argument_parser()
        args = parser.parse_args()
        # Figure out collection names
        if args.output_run is None:
            if args.output is None:
                raise ValueError("At least one of --output or --output-run options is required.")
            args.output_run = "{}/{}".format(args.output, Instrument.makeCollectionTimestamp())

        butler = Butler(args.repo, collections=args.collections)
        builder = cls(
            butler,
            dataset_type_names=args.dataset_types,
            where=args.where,
            input_collections=args.collections,
            output_run=args.output_run,
        )
        qg: QuantumGraph = builder.build(
            # Metadata includes a subset of attributes defined in CmdLineFwk.
            metadata={
                "input": args.collections,
                "butler_argument": args.repo,
                "output": args.output,
                "output_run": args.output_run,
                "data_query": args.where,
                "time": f"{datetime.datetime.now()}",
            }
        )
        qg.saveUri(args.filename)
