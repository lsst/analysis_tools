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

"""Base class implementation for the classes needed in creating `PipelineTasks`
which execute `AnalysisTools`.

The classes defined in this module have all the required behaviors for
defining, introspecting, and executing `AnalysisTools` against an input dataset
type.

Subclasses of these tasks should specify specific datasets to consume in their
connection classes and should specify a unique name
"""

__all__ = ("AnalysisBaseConfig", "AnalysisPipelineTask")

from collections import abc
from typing import TYPE_CHECKING, Any, Iterable, Mapping, MutableMapping, cast

if TYPE_CHECKING:
    from lsst.daf.butler import DeferredDatasetHandle

from lsst.daf.butler import DataCoordinate
from lsst.pex.config import ListField
from lsst.pex.config.configurableActions import ConfigurableActionStructField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base import connectionTypes as ct
from lsst.pipe.base.butlerQuantumContext import ButlerQuantumContext
from lsst.pipe.base.connections import InputQuantizedConnection, OutputQuantizedConnection

from ..analysisMetrics.metricMeasurementBundle import MetricMeasurementBundle
from ..interfaces import AnalysisMetric, AnalysisPlot, KeyedData


class AnalysisBaseConnections(
    PipelineTaskConnections, dimensions={}, defaultTemplates={"outputName": "Placeholder"}
):
    r"""Base class for Connections used for AnalysisTools PipelineTasks.

    This class has a pre-defined output connection for the
    MetricMeasurementMapping. The dataset type name for this connection is
    determined by the template ``outputName``.

    Output connections for plots created by `AnalysisPlot`\ s are created
    dynamically when an instance of the class is created. The init method
    examines all the `AnalysisPlot` actions specified in the associated
    `AnalysisBaseConfig` subclass accumulating all the info needed to
    create the output connections.

    The dimensions for all of the output connections (metric and plot) will
    be the same as the dimensions specified for the AnalysisBaseConnections
    subclass (i.e. quantum dimensions).
    """

    metrics = ct.Output(
        doc="Metrics calculated on input dataset type",
        name="{outputName}_metrics",
        storageClass="MetricMeasurementBundle",
    )

    def __init__(self, *, config: AnalysisBaseConfig = None):  # type: ignore
        # Validate that the outputName template has been set in config. This
        # should have been checked early with the configs validate method, but
        # it is possible for someone to manually create everything in a script
        # without running validate, so also check it late here.
        if (outputName := config.connections.outputName) == "Placeholder":  # type: ignore
            raise RuntimeError(
                "Subclasses must specify an alternative value for the defaultTemplate `outputName`"
            )
        super().__init__(config=config)

        # All arguments must be passed by kw, but python has not method to do
        # that without specifying a default, so None is used. Validate that
        # it is not None. This is largely for typing reasons, as in the normal
        # course of operation code execution paths ensure this will not be None
        assert config is not None

        # Set the dimensions for the metric
        self.metrics = ct.Output(
            name=self.metrics.name,
            doc=self.metrics.doc,
            storageClass=self.metrics.storageClass,
            dimensions=self.dimensions,
            multiple=False,
            isCalibration=False,
        )

        # Look for any conflicting names, creating a set of them, as these
        # will be added to the instance as well as recorded in the outputs
        # set.
        existingNames = set(dir(self))

        # Accumulate all the names to be used from all of the defined
        # AnalysisPlots.
        names: list[str] = []
        for plotAction in config.plots:
            if plotAction.parameterizedBand:
                for band in config.bands:
                    names.extend(name.format(band=band) for name in plotAction.getOutputNames())
            else:
                names.extend(plotAction.getOutputNames())

        # For each of the names found, create output connections.
        for name in names:
            name = f"{outputName}_{name}"
            if name in self.outputs or name in existingNames:
                raise NameError(
                    f"Plot with name {name} conflicts with existing connection"
                    " are two plots named the same?"
                )
            outConnection = ct.Output(
                name=name,
                storageClass="Plot",
                doc="Dynamic connection for plotting",
                dimensions=self.dimensions,
            )
            object.__setattr__(self, name, outConnection)
            self.outputs.add(name)


class AnalysisBaseConfig(PipelineTaskConfig, pipelineConnections=AnalysisBaseConnections):
    """Base class for all configs used to define an `AnalysisPipelineTask`

    This base class defines three fields that should be used in all subclasses,
    plots, metrics, and bands.

    The ``plots`` field is where a user configures which `AnalysisPlots` will
    be run in this `PipelineTask`.

    Likewise ``metrics`` defines which `AnalysisMetrics` will be run.

    The bands field specifies which bands will be looped over for
    `AnalysisTools` which support parameterized bands. I.e. called once for
    each band in the list.
    """

    plots = ConfigurableActionStructField[AnalysisPlot](doc="AnalysisPlots to run with this Task")
    metrics = ConfigurableActionStructField[AnalysisMetric](doc="AnalysisMetrics to run with this Task")
    bands = ListField[str](
        doc="Filter bands on which to run all of the actions", default=["u", "g", "r", "i", "z", "y"]
    )

    def validate(self):
        super().validate()
        # Validate that the required connections template is set.
        if self.connections.outputName == "Placeholder":  # type: ignore
            raise RuntimeError("Connections class 'outputName' must have a config explicitly set")


class _StandinPlotInfo(dict):
    """This class is an implementation detail to support plots in the instance
    no PlotInfo object is present in the call to run.
    """

    def __missing__(self, key):
        return ""


class AnalysisPipelineTask(PipelineTask):
    """Base class for `PipelineTasks` intended to run `AnalysisTools`.

    The run method will run all of the `AnalysisMetrics` and `AnalysisPlots`
    defined in the config class.

    To support interactive investigations, the actual work is done in
    ``runMetrics`` and ``runPlots`` methods. These can be called interactively
    with the same arguments as ``run`` but only the corresponding outputs will
    be produced.
    """

    # Typing config because type checkers dont know about our Task magic
    config: AnalysisBaseConfig
    ConfigClass = AnalysisBaseConfig

    def runPlots(self, data: KeyedData, **kwargs) -> Struct:
        results = Struct()
        # allow not sending in plot info
        if "plotInfo" not in kwargs:
            kwargs["plotInfo"] = _StandinPlotInfo()
        for name, action in self.config.plots.items():
            for selector in action.prep.selectors:
                if "threshold" in selector.keys():
                    kwargs["plotInfo"]["SN"] = selector.threshold
            kwargs["plotInfo"]["plotName"] = name
            match action(data, **kwargs):
                case abc.Mapping() as val:
                    for n, v in val.items():
                        setattr(results, n, v)
                case value:
                    setattr(results, name, value)
        if "SN" not in kwargs["plotInfo"].keys():
            kwargs["plotInfo"]["SN"] = "-"
        return results

    def runMetrics(self, data: KeyedData, **kwargs) -> Struct:
        metricsMapping = MetricMeasurementBundle()
        for name, action in self.config.metrics.items():
            match action(data, **kwargs):
                case abc.Mapping() as val:
                    results = list(val.values())
                case val:
                    results = [val]
            metricsMapping[name] = results  # type: ignore
        return Struct(metrics=metricsMapping)

    def run(self, *, data: KeyedData | None = None, **kwargs) -> Struct:
        """Produce the outputs associated with this `PipelineTask`

        Parameters
        ----------
        data : `KeyedData`
            The input data from which all `AnalysisTools` will run and produce
            outputs. A side note, the python typing specifies that this can be
            None, but this is only due to a limitation in python where in order
            to specify that all arguments be passed only as keywords the
            argument must be given a default. This argument most not actually
            be None.
        **kwargs
            Additional arguments that are passed through to the `AnalysisTools`
            specified in the configuration.

        Returns
        -------
        results : `~lsst.pipe.base.Struct`
            The accumulated results of all the plots and metrics produced by
            this `PipelineTask`.

        Raises
        ------
        ValueError
            Raised if the supplied data argument is `None`
        """
        if data is None:
            raise ValueError("data must not be none")
        results = Struct()
        plotKey = f"{self.config.connections.outputName}_{{name}}"  # type: ignore
        if "bands" not in kwargs:
            kwargs["bands"] = list(self.config.bands)
        kwargs["plotInfo"]["bands"] = kwargs["bands"]
        for name, value in self.runPlots(data, **kwargs).getDict().items():
            setattr(results, plotKey.format(name=name), value)
        for name, value in self.runMetrics(data, **kwargs).getDict().items():
            setattr(results, name, value)

        return results

    def runQuantum(
        self,
        butlerQC: ButlerQuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        """Override default runQuantum to load the minimal columns necessary
        to complete the action.

        Parameters
        ----------
        butlerQC : `ButlerQuantumContext`
            A butler which is specialized to operate in the context of a
            `lsst.daf.butler.Quantum`.
        inputRefs : `InputQuantizedConnection`
            Datastructure whose attribute names are the names that identify
            connections defined in corresponding `PipelineTaskConnections`
            class. The values of these attributes are the
            `lsst.daf.butler.DatasetRef` objects associated with the defined
            input/prerequisite connections.
        outputRefs : `OutputQuantizedConnection`
            Datastructure whose attribute names are the names that identify
            connections defined in corresponding `PipelineTaskConnections`
            class. The values of these attributes are the
            `lsst.daf.butler.DatasetRef` objects associated with the defined
            output connections.
        """
        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId
        plotInfo = self.parsePlotInfo(inputs, dataId)
        data = self.loadData(inputs["data"])
        if "skymap" in inputs.keys():
            skymap = inputs["skymap"]
        else:
            skymap = None
        outputs = self.run(data=data, plotInfo=plotInfo, skymap=skymap)
        butlerQC.put(outputs, outputRefs)

    def _populatePlotInfoWithDataId(
        self, plotInfo: MutableMapping[str, Any], dataId: DataCoordinate | None
    ) -> None:
        """Update the plotInfo with the dataId values.

        Parameters
        ----------
        plotInfo : `dict`
            The plotInfo dictionary to update.
        dataId : `lsst.daf.butler.DataCoordinate`
            The dataId to use to update the plotInfo.
        """
        if dataId is not None:
            for dataInfo in dataId:
                plotInfo[dataInfo.name] = dataId[dataInfo.name]

    def parsePlotInfo(
        self, inputs: Mapping[str, Any] | None, dataId: DataCoordinate | None, connectionName: str = "data"
    ) -> Mapping[str, str]:
        """Parse the inputs and dataId to get the information needed to
        to add to the figure.

        Parameters
        ----------
        inputs: `dict`
            The inputs to the task
        dataCoordinate: `lsst.daf.butler.DataCoordinate`
            The dataId that the task is being run on.
        connectionName: `str`, optional
            Name of the input connection to use for determining table name.

        Returns
        -------
        plotInfo : `dict`
        """

        if inputs is None:
            tableName = ""
            run = ""
        else:
            tableName = inputs[connectionName].ref.datasetType.name
            run = inputs[connectionName].ref.run

        # Initialize the plot info dictionary
        plotInfo = {"tableName": tableName, "run": run}

        self._populatePlotInfoWithDataId(plotInfo, dataId)
        return plotInfo

    def loadData(self, handle: DeferredDatasetHandle, names: Iterable[str] | None = None) -> KeyedData:
        """Load the minimal set of keyed data from the input dataset.

        Parameters
        ----------
        handle : `DeferredDatasetHandle`
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
        if names is None:
            names = self.collectInputNames()
        return cast(KeyedData, handle.get(parameters={"columns": names}))

    def collectInputNames(self) -> Iterable[str]:
        """Get the names of the inputs.

        If using the default `loadData` method this will gather the names
        of the keys to be loaded from an input dataset.

        Returns
        -------
        inputs : `Iterable` of `str`
            The names of the keys in the `KeyedData` object to extract.

        """
        inputs = set()
        for band in self.config.bands:
            for name, action in self.config.plots.items():
                for column, dataType in action.getFormattedInputSchema(band=band):
                    inputs.add(column)
            for name, action in self.config.metrics.items():
                for column, dataType in action.getFormattedInputSchema(band=band):
                    inputs.add(column)
        return inputs
