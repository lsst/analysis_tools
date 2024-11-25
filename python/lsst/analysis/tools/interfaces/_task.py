# This file is part of analysis_tools.  #
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

__all__ = ("AnalysisBaseConnections", "AnalysisBaseConfig", "AnalysisPipelineTask")

import datetime
import logging
import warnings
import weakref
from collections.abc import Iterable
from copy import deepcopy
from typing import TYPE_CHECKING, Any, Mapping, MutableMapping, cast

import matplotlib.pyplot as plt
from lsst.verify import Measurement

if TYPE_CHECKING:
    from lsst.daf.butler import DeferredDatasetHandle
    from lsst.pipe.base import QuantumContext

from lsst.daf.butler import DataCoordinate
from lsst.pex.config import Field, ListField
from lsst.pex.config.configurableActions import ConfigurableActionStructField
from lsst.pipe.base import Instrument, PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base import connectionTypes as ct
from lsst.pipe.base.connections import InputQuantizedConnection, OutputQuantizedConnection
from lsst.pipe.base.pipelineIR import ConfigIR, ParametersIR

from ._actions import JointAction, MetricAction, NoMetric
from ._analysisTools import AnalysisTool
from ._interfaces import KeyedData, PlotTypes
from ._metricMeasurementBundle import MetricMeasurementBundle

# TODO: This rcParams modification is a temporary solution, hiding
# a matplotlib warning indicating too many figures have been opened.
# When DM-39114 is implemented, this should be removed.
plt.rcParams.update({"figure.max_open_warning": 0})


# TODO: This _plotCloser function assists in closing all open plots at the
# conclusion of a PipelineTask. When DM-39114 is implemented, this function and
# all associated usage thereof should be removed.
def _plotCloser(*args):
    """Close all the plots in the given list."""
    for plot in args:
        plt.close(plot)


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

        if config.dimensions:
            self.dimensions = set(config.dimensions)

        for tool in config.atools:
            match tool.produce:
                case JointAction():
                    if isinstance(tool.produce.metric, NoMetric):
                        continue
                    if len(tool.produce.metric.units) != 0:
                        hasMetrics = True
                        break
                case MetricAction():
                    hasMetrics = True
                    break
        else:
            hasMetrics = False

        # Set the dimensions for the metric
        if hasMetrics and False:
            self.metrics = ct.Output(
                name=self.metrics.name,
                doc=self.metrics.doc,
                storageClass=self.metrics.storageClass,
                dimensions=self.dimensions,
                multiple=False,
                isCalibration=False,
            )
        else:
            # There are no metrics to produce, remove the output connection
            self.outputs.remove("metrics")

        # Look for any conflicting names, creating a set of them, as these
        # will be added to the instance as well as recorded in the outputs
        # set.
        existingNames = set(dir(self))

        # Accumulate all the names to be used from all of the defined
        # AnalysisPlots.
        # names: Mapping[str, ConfigurableAction] = {}
        names = {}
        for action in config.atools:
            if action.dynamicOutputNames:
                outNames = action.getOutputNames(config=config)
            else:
                outNames = action.getOutputNames()
            if action.parameterizedBand and not "band" in self.dimensions:
                for band in config.bands:
                    # names.extend(name.format(band=band) for name in outNames)
                    # names.extend("_".join((band, name)) for name in outNames)
                    names.update({name.format(band=band): action for name in outNames})
            else:
                # names.extend(outNames)
                names.update({name: action for name in outNames})

        # For each of the names found, create output connections.
        for name, action in names.items():
            name = f"{outputName}_{name}"
            if name in self.outputs or name in existingNames:
                raise NameError(
                    f"Plot with name {name} conflicts with existing connection"
                    " are two plots named the same?"
                )

            if action.parameterizedBand and "band" not in self.dimensions:
                multiple = True
                dimensions = self.dimensions.union({"band"})
            else:
                multiple = False
                dimensions = self.dimensions

            outConnection = ct.Output(
                name=name,
                storageClass="Plot",
                doc="Dynamic connection for plotting",
                dimensions=dimensions,
                multiple=multiple,
            )
            setattr(self, name, outConnection)


def _timestampValidator(value: str) -> bool:
    if value in ("reference_package_timestamp", "run_timestamp", "current_timestamp", "dataset_timestamp"):
        return True
    elif "explicit_timestamp" in value:
        try:
            _, splitTime = value.split(":")
        except ValueError:
            logging.error(
                "Explicit timestamp must be given in the format 'explicit_timestamp:datetime', "
                r"where datetime is given in the form '%Y%m%dT%H%M%S%z"
            )
            return False
        try:
            datetime.datetime.strptime(splitTime, r"%Y%m%dT%H%M%S%z")
        except ValueError:
            # This is explicitly chosen to be an f string as the string
            # contains control characters.
            logging.error(
                f"The supplied datetime {splitTime} could not be parsed correctly into "
                r"%Y%m%dT%H%M%S%z format"
            )
            return False
        return True
    else:
        return False


class AnalysisBaseConfig(PipelineTaskConfig, pipelineConnections=AnalysisBaseConnections):
    """Base class for all configs used to define an `AnalysisPipelineTask`.

    This base class defines two fields that should be used in all subclasses,
    atools, and bands.

    The ``atools`` field is where the user configures which analysis tools will
    be run as part of this `PipelineTask`.

    The bands field specifies which bands will be looped over for
    `AnalysisTools` which support parameterized bands. I.e. called once for
    each band in the list.
    """

    atools = ConfigurableActionStructField[AnalysisTool](
        doc="The analysis tools that are to be run by this task at execution"
    )
    # Temporarally alias these for backwards compatibility
    plots = atools
    metrics = atools
    bands = ListField[str](
        doc="Filter bands on which to run all of the actions",
        default=[],
        optional=True,
    )
    metric_tags = ListField[str](
        doc="List of tags which will be added to all configurable actions", default=[]
    )
    dataset_identifier = Field[str](doc="An identifier to be associated with output Metrics", optional=True)
    reference_package = Field[str](
        doc="A package who's version, at the time of metric upload to a "
        "time series database, will be converted to a timestamp of when "
        "that version was produced",
        default="lsst_distrib",
    )
    timestamp_version = Field[str](
        doc="Which time stamp should be used as the reference timestamp for a "
        "metric in a time series database, valid values are; "
        "reference_package_timestamp, run_timestamp, current_timestamp, "
        "dataset_timestamp and explicit_timestamp:datetime where datetime is "
        "given in the form %Y%m%dT%H%M%S%z",
        default="run_timestamp",
        check=_timestampValidator,
    )
    dimensions = ListField[str](
        doc="The dimensions of the input and output connections", default=None, optional=True
    )

    def applyConfigOverrides(
        self,
        instrument: Instrument | None,
        taskDefaultName: str,
        pipelineConfigs: Iterable[ConfigIR] | None,
        parameters: ParametersIR,
        label: str,
    ) -> None:
        extraConfig = {}
        if (value := parameters.mapping.get("sasquatch_dataset_identifier", None)) is not None:
            extraConfig["dataset_identifier"] = value
        if (value := parameters.mapping.get("sasquatch_reference_package", None)) is not None:
            extraConfig["reference_package"] = value
        if (value := parameters.mapping.get("sasquatch_timestamp_version", None)) is not None:
            if "explicit_timestamp" in value:
                try:
                    _, splitTime = value.split(":")
                except ValueError as excpt:
                    raise ValueError(
                        "Explicit timestamp must be given in the format 'explicit_timestamp:datetime', "
                        "where datetime is given in the form '%Y%m%dT%H%M%S%z"
                    ) from excpt
                try:
                    datetime.datetime.strptime(splitTime, r"%Y%m%dT%H%M%S%z")
                except ValueError as excpt:
                    raise ValueError(
                        f"The supplied datetime {splitTime} could not be parsed correctly into "
                        "%Y%m%dT%H%M%S%z format"
                    ) from excpt
            extraConfig["timestamp_version"] = value
        if extraConfig:
            newPipelineConfigs = [ConfigIR(rest=extraConfig)]
            if pipelineConfigs is not None:
                newPipelineConfigs.extend(pipelineConfigs)
            pipelineConfigs = newPipelineConfigs
        return super().applyConfigOverrides(instrument, taskDefaultName, pipelineConfigs, parameters, label)

    def freeze(self):
        # Copy the meta configuration values to each of the configured tools
        # only do this if the tool has not been further specialized
        if not self._frozen:
            for tool in self.atools:
                for tag in self.metric_tags:
                    tool.metric_tags.insert(-1, tag)
        super().freeze()

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

    The run method will run all of the `AnalysisTools` defined in the config
    class.
    """

    # Typing config because type checkers dont know about our Task magic
    config: AnalysisBaseConfig
    ConfigClass = AnalysisBaseConfig

    warnings_all = (
        "divide by zero encountered in divide",
        "invalid value encountered in arcsin",
        "invalid value encountered in cos",
        "invalid value encountered in divide",
        "invalid value encountered in log10",
        "invalid value encountered in scalar divide",
        "invalid value encountered in sin",
        "invalid value encountered in sqrt",
        "invalid value encountered in true_divide",
        "Mean of empty slice",
    )

    def _runTools(self, data: KeyedData, **kwargs) -> Struct:
        with warnings.catch_warnings():
            # Change below to "in self.warnings_all" to find otherwise
            # unfiltered numpy warnings.
            for warning in ():
                warnings.filterwarnings("error", warning, RuntimeWarning)
            results = Struct()
            results.metrics = MetricMeasurementBundle(
                dataset_identifier=self.config.dataset_identifier,
                reference_package=self.config.reference_package,
                timestamp_version=self.config.timestamp_version,
            )
            # copy plot info to be sure each action sees its own copy
            plotInfo = kwargs.get("plotInfo")
            plotKey = f"{self.config.connections.outputName}_{{name}}"
            weakrefArgs = []
            for name, action in self.config.atools.items():
                kwargs["plotInfo"] = deepcopy(plotInfo)
                actionResult = action(data, **kwargs)
                metricAccumulate = []
                for resultName, value in actionResult.items():
                    match value:
                        case PlotTypes():
                            setattr(results, plotKey.format(name=resultName), value)
                            weakrefArgs.append(value)
                        case Measurement():
                            metricAccumulate.append(value)
                # only add the metrics if there are some
                if metricAccumulate:
                    results.metrics[name] = metricAccumulate
            # Wrap the return struct in a finalizer so that when results is
            # garbage collected the plots will be closed.
            # TODO: This finalize step closes all open plots at the conclusion
            # of a task. When DM-39114 is implemented, this step should not
            # be required and may be removed.
            weakref.finalize(results, _plotCloser, *weakrefArgs)
        return results

    def run(self, *, data: KeyedData | None = None, **kwargs) -> Struct:
        """Produce the outputs associated with this `PipelineTask`.

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
        if "bands" not in kwargs:
            kwargs["bands"] = list(self.config.bands)
        if "plotInfo" not in kwargs:
            kwargs["plotInfo"] = _StandinPlotInfo()
        kwargs["plotInfo"]["bands"] = kwargs["bands"]
        return self._runTools(data, **kwargs)

    def runQuantum(
        self,
        butlerQC: QuantumContext,
        inputRefs: InputQuantizedConnection,
        outputRefs: OutputQuantizedConnection,
    ) -> None:
        """Override default runQuantum to load the minimal columns necessary
        to complete the action.

        Parameters
        ----------
        butlerQC : `~lsst.pipe.base.QuantumContext`
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
        # We implicitly assume that 'data' has been defined, but do not have a
        # corresponding input connection in the base class. Thus, we capture
        # and re-raise the error with a more helpful message.
        try:
            # data has to be popped out to avoid duplication in the call to the
            # `run` method.
            inputData = inputs.pop("data")
        except KeyError:
            raise RuntimeError("'data' is a required input connection, but is not defined.")
        if band := dataId.get("band", None):
            inputs["band"] = band
        data = self.loadData(inputData, dataId=dataId)
        outputs = self.run(data=data, plotInfo=plotInfo, **inputs)

        if "band" in dataId:
            # Branch for when the PipelineTask has band in its dimensions.
            butlerQC.put(outputs, outputRefs)
        else:
            for outputRefName in outputRefs.keys():
                name = outputRefName.split(self.config.connections.outputName, maxsplit=1)[1]
                # _, name = band_name.split("_", maxsplit=2)[1:]
                if outputRefName == "metrics":
                    butlerQC.put(getattr(outputs, outputRefName), getattr(outputRefs, outputRefName))
                    continue

                try:
                    datasetRef = getattr(outputRefs, outputRefName)
                    if hasattr(datasetRef, "__iter__"):
                        for outputRef in datasetRef:
                            band = outputRef.dataId["band"]
                            # name would already have a leading underscore.
                            newOutputName = f"{self.config.connections.outputName}_{band}{name}"
                            if dataset := getattr(outputs, newOutputName, None):
                                butlerQC.put(dataset, outputRef)
                    else:
                        if dataset := getattr(outputs, outputRefName, None):
                            butlerQC.put(dataset, datasetRef)
                except Exception as e:
                    pass
                    raise e

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
            plotInfo.update(dataId.mapping)

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

    def loadData(
        self, handle: DeferredDatasetHandle, names: Iterable[str] | None = None, dataId=None
    ) -> KeyedData:
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
            names = self.collectInputNames(dataId=dataId)
        return cast(KeyedData, handle.get(parameters={"columns": names}))

    def collectInputNames(self, dataId=None) -> Iterable[str]:
        """Get the names of the inputs.

        If using the default `loadData` method this will gather the names
        of the keys to be loaded from an input dataset.

        Returns
        -------
        inputs : `Iterable` of `str`
            The names of the keys in the `KeyedData` object to extract.

        """
        inputs = set()

        if not (localBands := self.config.bands):
            localBands = set("")
        if dataId is not None:
            localBands.update(dataId.get("band", ""))
        for band in localBands:
            for action in self.config.atools:
                for key, _ in action.getFormattedInputSchema(band=band):
                    inputs.add(key)
        return inputs
