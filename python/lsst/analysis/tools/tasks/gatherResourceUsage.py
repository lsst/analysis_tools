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
)

import argparse
import dataclasses
import datetime
import logging
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
import pandas as pd
from lsst.daf.butler import Butler, DataCoordinate, DatasetRef, DatasetType, DimensionGraph, Quantum
from lsst.daf.butler.core.utils import globToRegex
from lsst.pex.config import Field, ListField
from lsst.pipe.base import (
    Instrument,
    PipelineDatasetTypes,
    PipelineTask,
    PipelineTaskConfig,
    PipelineTaskConnections,
    QuantumGraph,
    Struct,
    TaskDef,
)
from lsst.pipe.base import connectionTypes as cT
from lsst.utils.introspection import get_full_type_name

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
        memrun = pd.merge(
            memoryGB.reset_index(),
            runtime.reset_index(),
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
        doc=("Whether to extract the CPU time duration for actually " "constructing the task."),
        default=True,
    )
    run_time = Field[bool](
        doc=("Whether to extract the CPU time duration for actually " "executing the task."),
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

    @classmethod
    def build_quantum_graph(
        cls, metadata_refs: Iterable[DatasetRef], graph_metadata: Mapping[str, Any]
    ) -> QuantumGraph:
        """Build a specialized `QuantumGraph` that configures and runs this
        task on existing metadata datasets.

        Parameters
        ----------
        metadata_refs : `Iterable` [ `DatasetRef` ]
            References to metadata datasets.  Non-metadata datasets are
            silently ignored, in order to make it easier to pass in more
            general dataset query results.
        graph_metadata : `dict` [ `str`, `Any` ]
            Graph metadata. It is required to contain "output_run" key with the
            name of the output RUN collection.

        Returns
        -------
        qg : `lsst.pipe.base.QuantumGraph`
            New quantum graph.

        Notes
        -----
        This task cannot easily be added to a regular pipeline, as it's much
        more natural to have it run automatically on all *other* tasks.
        But even machine-generating a pipeline is problematic, because our
        current `QuantumGraph` generation algorithm isn't smart enough to
        recognize that those pipelines are usually disconnected subgraphs, and
        that generating a graph for each of those separately is many orders of
        magnitude faster (and *not* doing is disastrously slow).  As a
        short-term workaround, this method can be used to generate a
        `QuantumGraph` for just this task much more efficiently from the
        metadata datasets already present in a set of collections.
        """
        # Group input metadata datasets by dataset type.
        metadata_refs_by_dataset_type: dict[DatasetType, set[DatasetRef]] = defaultdict(set)
        for ref in metadata_refs:
            metadata_refs_by_dataset_type[ref.datasetType].add(ref)
        metadata_refs_by_dataset_type = dict(metadata_refs_by_dataset_type)
        # Iterate over those groups, creating a configuration of
        # this task and quanta for each one.
        quanta_by_task_def = {}
        init_outputs = {}
        empty_dimensions: DimensionGraph | None = None
        data_id: DataCoordinate | None = None

        consolidate_inputs = {}
        consolidate_config = ConsolidateResourceUsageConfig()

        for input_metadata_dataset_type, metadata_refs in metadata_refs_by_dataset_type.items():
            if (m := re.fullmatch(r"^(\w+)_metadata$", input_metadata_dataset_type.name)) is None:
                continue
            elif "gatherResourceUsage" in input_metadata_dataset_type.name:
                continue
            else:
                input_task_label = m.group(1)
            _LOG.info(
                "Creating GatherResourceUsage quantum for %s with %d input datasets.",
                input_task_label,
                len(metadata_refs),
            )
            dataset_type_name = f"{input_task_label}_resource_usage"

            config = cls.ConfigClass()
            config.dimensions = input_metadata_dataset_type.dimensions.names
            config.connections.input_metadata = input_metadata_dataset_type.name
            config.connections.output_table = dataset_type_name

            consolidate_config.input_names.append(dataset_type_name)

            task_def = TaskDef(
                taskName=get_full_type_name(cls),
                taskClass=cls,
                config=config,
                label=f"{input_task_label}_gatherResourceUsage",
            )
            empty_dimensions = input_metadata_dataset_type.dimensions.universe.empty
            output_table_dataset_type = DatasetType(
                config.connections.output_table,
                dimensions=empty_dimensions,
                storageClass=GatherResourceUsageConnections.output_table.storageClass,
            )
            data_id = DataCoordinate.makeEmpty(universe=input_metadata_dataset_type.dimensions.universe)
            output_run = graph_metadata["output_run"]
            outputs = {
                output_table_dataset_type: [DatasetRef(output_table_dataset_type, data_id, run=output_run)],
            }
            consolidate_inputs.update(outputs)
            if task_def.metadataDatasetName is not None:
                output_metadata_dataset_type = DatasetType(
                    task_def.metadataDatasetName,
                    dimensions=empty_dimensions,
                    storageClass="TaskMetadata",
                )
                outputs[output_metadata_dataset_type] = [
                    DatasetRef(output_metadata_dataset_type, data_id, run=output_run)
                ]
            if task_def.logOutputDatasetName is not None:
                log_dataset_type = DatasetType(
                    task_def.logOutputDatasetName,
                    dimensions=empty_dimensions,
                    storageClass="ButlerLogRecords",
                )
                outputs[log_dataset_type] = [DatasetRef(log_dataset_type, data_id, run=output_run)]
            quantum = Quantum(
                taskName=task_def.taskName,
                taskClass=task_def.taskClass,
                dataId=data_id,
                inputs={input_metadata_dataset_type: list(metadata_refs)},
                outputs=outputs,
            )
            quanta_by_task_def[task_def] = {quantum}

            config_dataset_type = DatasetType(
                task_def.configDatasetName,
                dimensions=empty_dimensions,
                storageClass="Config",
            )
            init_outputs[task_def] = [DatasetRef(config_dataset_type, data_id, run=output_run)]

        if empty_dimensions is None:
            raise RuntimeError("No metadata dataset refs given.")
        assert data_id is not None

        # Now once more for the task to consolidate all the individual resource
        # usage tables
        task_def = TaskDef(
            taskName=get_full_type_name(ConsolidateResourceUsageTask),
            taskClass=ConsolidateResourceUsageTask,
            config=consolidate_config,
            label=ConsolidateResourceUsageTask._DefaultName,
        )

        output_table_dataset_type = DatasetType(
            consolidate_config.connections.output_table,
            dimensions=empty_dimensions,
            storageClass=ConsolidateResourceUsageConnections.output_table.storageClass,
        )

        outputs = {
            output_table_dataset_type: [DatasetRef(output_table_dataset_type, data_id, run=output_run)],
        }

        if task_def.metadataDatasetName is not None:
            output_metadata_dataset_type = DatasetType(
                task_def.metadataDatasetName,
                dimensions=empty_dimensions,
                storageClass="TaskMetadata",
            )
            outputs[output_metadata_dataset_type] = [
                DatasetRef(output_metadata_dataset_type, data_id, run=output_run)
            ]

        if task_def.logOutputDatasetName is not None:
            log_dataset_type = DatasetType(
                task_def.logOutputDatasetName,
                dimensions=empty_dimensions,
                storageClass="ButlerLogRecords",
            )
            outputs[log_dataset_type] = [DatasetRef(log_dataset_type, data_id, run=output_run)]

        quantum = Quantum(
            taskName=task_def.taskName,
            taskClass=task_def.taskClass,
            dataId=data_id,
            inputs=consolidate_inputs,
            outputs=outputs,
        )
        quanta_by_task_def[task_def] = {quantum}

        global_init_outputs = []
        if empty_dimensions is not None and data_id is not None:
            packages_dataset_type = DatasetType(
                PipelineDatasetTypes.packagesDatasetName,
                dimensions=empty_dimensions,
                storageClass="Packages",
            )
            global_init_outputs.append(DatasetRef(packages_dataset_type, data_id, run=output_run))

        return QuantumGraph(
            quanta_by_task_def,
            initOutputs=init_outputs,
            globalInitOutputs=global_init_outputs,
            metadata=graph_metadata,
        )

    @classmethod
    def build_quantum_graph_cli(cls, argv):
        """Run the command-line interface for `build_quantum_graph`.

        This method provides the implementation for the
        ``build-gather-resource-usage-qg`` script.

        Parameters
        ----------
        argv : `Sequence` [ `str` ]
            Command-line arguments (e.g. ``sys.argv[1:]``).
        """
        parser = argparse.ArgumentParser(
            description=(
                "Build a QuantumGraph that runs GatherResourceUsageTask on existing metadata datasets."
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
            help=(
                "Glob-style patterns for input metadata dataset types.  If a pattern matches a "
                "non-metadata dataset type, non-metadata matches will be ignored."
            ),
        )
        parser.add_argument(
            "--where",
            type=str,
            default=None,
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
        args = parser.parse_args(argv)
        base_dataset_type_filter = re.compile(r"\w+_metadata")
        if not args.dataset_types:
            input_dataset_types = base_dataset_type_filter
        else:
            input_dataset_types = [globToRegex(expr) for expr in args.dataset_types]
        butler = Butler(args.repo, collections=args.collections)
        metadata_refs = butler.registry.queryDatasets(input_dataset_types, where=args.where, findFirst=True)

        # Figure out collection names
        if args.output_run is None:
            if args.output is None:
                raise ValueError("At least one of --output or --output-run options is required.")
            args.output_run = "{}/{}".format(args.output, Instrument.makeCollectionTimestamp())

        # Metadata includes a subset of attributes defined in CmdLineFwk.
        graph_metadata = {
            "input": args.collections,
            "butler_argument": args.repo,
            "output": args.output,
            "output_run": args.output_run,
            "data_query": args.where,
            "time": f"{datetime.datetime.now()}",
        }
        qg = cls.build_quantum_graph(metadata_refs, graph_metadata)
        qg.saveUri(args.filename)

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
        dimensions = universe.extract(self.config.dimensions)
        # Transform input list into a dict keyed by data ID.
        handles_by_data_id = {}
        for handle in input_metadata:
            handles_by_data_id[handle.dataId] = handle
        n_rows = len(handles_by_data_id)
        # Create a dict of empty column arrays that we'll ultimately make into
        # a table.
        columns = {d.name: np.zeros(n_rows, dtype=_dtype_from_field_spec(d.primaryKey)) for d in dimensions}
        for attr_name in ("memory", "prep_time", "init_time", "run_time"):
            if getattr(self.config, attr_name):
                columns[attr_name] = np.zeros(n_rows, dtype=float)
        for method_name in self.config.method_times:
            columns[method_name] = np.zeros(n_rows, dtype=float)
        # Populate the table, one row at a time.
        warned_about_metadata_version = False
        for index, (data_id, handle) in enumerate(handles_by_data_id.items()):
            # Fill in the data ID columns.
            for k, v in data_id.full.byName().items():
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
        return {
            attr_name: end - begin
            for attr_name, begin, end in zip(
                ["prep_time", "init_time", "run_time"],
                times[:-1],
                times[1:],
            )
            if getattr(self.config, attr_name)
        }

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
