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

__all__ = ("AnalysisTool",)

from collections import ChainMap
from collections.abc import Mapping
from functools import wraps
from operator import attrgetter
from typing import Callable, Iterable, Protocol, runtime_checkable

import lsst.pex.config as pexConfig
from lsst.obs.base import Instrument
from lsst.pex.config import Field, FieldValidationError, ListField
from lsst.pex.config.configurableActions import ConfigurableActionField
from lsst.pipe.base import Pipeline
from lsst.verify import Measurement

from ._actions import AnalysisAction, JointAction, JointResults, NoPlot, PlotAction
from ._interfaces import KeyedData, KeyedDataSchema, KeyedResults, PlotTypes
from ._stages import BasePrep, BaseProcess, BaseProduce


@runtime_checkable
class _HasOutputNames(Protocol):
    def getOutputNames(self, config: pexConfig.Config | None = None) -> Iterable[str]: ...


def _finalizeWrapper(
    f: Callable[[AnalysisTool], None], cls: type[AnalysisTool]
) -> Callable[[AnalysisTool], None]:
    """Wrap a classes finalize function to ensure the base classes special
    finalize method only fires after the most derived finalize method.

    Parameters
    ----------
    f : `Callable`
        Function that is being wrapped
    cls : `type` of `AnalysisTool`
        The class which is having its function wrapped

    Returns
    -------
    function : `Callable`
        The new function which wraps the old
    """

    @wraps(f)
    def wrapper(self: AnalysisTool) -> None:
        # call the wrapped finalize function
        f(self)
        # get the method resolution order for the self variable
        mro = self.__class__.mro()

        # Find which class in the mro that last defines a finalize method
        # note that this is in the reverse order from the mro, because the
        # last class in an inheritance stack is the first in the mro (aka you
        # walk from the furthest child first.
        #
        # Also note that the most derived finalize method need not be the same
        # as the type of self, as that might inherit from a parent somewhere
        # between it and the furthest parent.
        mostDerived: type | None = None
        for klass in mro:
            # inspect the classes dictionary to see if it specifically defines
            # finalize. This is needed because normal lookup will go through
            # the mro, but this needs to be restricted to each class.
            if "finalize" in vars(klass):
                mostDerived = klass
                break

        # Find what stage in the MRO walking process the recursive function
        # call is in.
        this = super(cls, self).__thisclass__

        # If the current place in the MRO walking is also the class that
        # defines the most derived instance of finalize, then call the base
        # classes private finalize that must be called after everything else.
        if mostDerived is not None and this == mostDerived:
            self._baseFinalize()

    return wrapper


class AnalysisTool(AnalysisAction):
    r"""A tool which which calculates a single type of analysis on input data,
    though it may return more than one result.

    Although `AnalysisTool`\ s are considered a single type of analysis, the
    classes themselves can be thought of as a container. `AnalysisTool`\ s
    are aggregations of `AnalysisAction`\ s to form prep, process, and
    produce stages. These stages allow better reuse of individual
    `AnalysisActions` and easier introspection in contexts such as a notebook
    or interpreter.

    An `AnalysisTool` can be thought of an an individual configuration that
    specifies which `AnalysisAction` should run for each stage.

    The stages themselves are also configurable, allowing control over various
    aspects of the individual `AnalysisAction`\ s.
    """

    prep = ConfigurableActionField[AnalysisAction](doc="Action to run to prepare inputs", default=BasePrep)
    process = ConfigurableActionField[AnalysisAction](
        doc="Action to process data into intended form", default=BaseProcess
    )
    produce = ConfigurableActionField[AnalysisAction](
        doc="Action to perform any finalization steps", default=BaseProduce
    )
    metric_tags = ListField[str](
        doc="List of tags which will be associated with metric measurement(s)", default=[]
    )

    def __init_subclass__(cls: type[AnalysisTool], **kwargs):
        super().__init_subclass__(**kwargs)
        # Wrap all definitions of the finalize method in a special wrapper that
        # ensures that the bases classes private finalize is called last.
        if "finalize" in vars(cls):
            cls.finalize = _finalizeWrapper(cls.finalize, cls)

    dynamicOutputNames: bool | Field[bool] = False
    """Determines whether to grant the ``getOutputNames`` method access to
    config parameters.
    """

    parameterizedBand: bool | Field[bool] = True
    """Specifies if an `AnalysisTool` may parameterize a band within any field
    in any stage, or if the set of bands is already uniquely determined though
    configuration. I.e. can this `AnalysisTool` be automatically looped over to
    produce a result for multiple bands.
    """

    propagateData: bool | Field[bool] = False
    """If this value is set to True, the input data `KeyedData` will be passed
    to each stage in addition to the ``prep`` stage. Any keys created in a
    stage with the same key that exists in the input ``data`` will shadow that
    key/value.
    """

    def __call__(self, data: KeyedData, **kwargs) -> KeyedResults:
        bands = kwargs.pop("bands", None)
        if "plotInfo" in kwargs and kwargs.get("plotInfo") is not None:
            if "plotName" not in kwargs["plotInfo"] or kwargs["plotInfo"]["plotName"] is None:
                kwargs["plotInfo"]["plotName"] = self.identity
        if not self.parameterizedBand or bands is None:
            if "band" not in kwargs:
                # Some tasks require a "band" key for naming. This shouldn't
                # affect the results. DM-35813 should make this unnecessary.
                kwargs["band"] = "analysisTools"
            return self._call_single(data, **kwargs)
        results: KeyedResults = {}
        for band in bands:
            kwargs["band"] = band
            if "plotInfo" in kwargs:
                kwargs["plotInfo"]["bands"] = band
            subResult = self._call_single(data, **kwargs)
            for key, value in subResult.items():
                match value:
                    case PlotTypes():
                        results[f"{band}_{key}"] = value
                    case Measurement():
                        results[key] = value
        return results

    def _call_single(self, data: KeyedData, **kwargs) -> KeyedResults:
        # create a shallow copy of kwargs
        kwargs = dict(**kwargs)
        kwargs["metric_tags"] = list(self.metric_tags or ())
        prepped: KeyedData = self.prep(data, **kwargs)  # type: ignore
        if self.propagateData:
            prepped = ChainMap(data, prepped)
        processed: KeyedData = self.process(prepped, **kwargs)  # type: ignore
        if self.propagateData:
            processed = ChainMap(data, processed)
        finalized: (
            Mapping[str, PlotTypes] | PlotTypes | Mapping[str, Measurement] | Measurement | JointResults
        ) = self.produce(
            processed, **kwargs
        )  # type: ignore
        return self._process_single_results(finalized)

    def _getPlotType(self) -> str:
        match self.produce:
            case PlotAction():
                return self.produce.getPlotType()
            case JointAction(plot=NoPlot()):
                pass
            case JointAction(plot=plotter):
                return plotter.getPlotType()

        return ""

    def _process_single_results(
        self,
        results: Mapping[str, PlotTypes] | PlotTypes | Mapping[str, Measurement] | Measurement | JointResults,
    ) -> KeyedResults:
        accumulation = {}
        predicate = f"{self.identity}" if self.identity else ""
        match results:
            case Mapping():
                suffix = self._getPlotType()
                for key, value in results.items():
                    match value:
                        case Measurement():
                            iterable = (predicate, key)
                        case PlotTypes():
                            iterable = (predicate, key, suffix)
                        case _:
                            raise RuntimeError(f"Unexpected {key=}, {value=} from:\n{self=}")
                    refKey = "_".join(x for x in iterable if x)
                    accumulation[refKey] = value
            case PlotTypes():
                suffix = self._getPlotType()
                refKey = "_".join(x for x in (predicate, suffix) if x)
                accumulation[refKey] = results
            case Measurement():
                accumulation[f"{predicate}"] = results
            case JointResults(plot=plotResults, metric=metricResults):
                if plotResults is not None:
                    subResult = self._process_single_results(plotResults)
                    accumulation.update(subResult)
                if metricResults is not None:
                    subResult = self._process_single_results(metricResults)
                    accumulation.update(subResult)
        return accumulation

    def getInputSchema(self) -> KeyedDataSchema:
        return self.prep.getInputSchema()

    def populatePrepFromProcess(self):
        """Add additional inputs to the prep stage if supported.

        If the configured prep action supports adding to its input schema,
        attempt to add the required inputs schema from the process stage to the
        prep stage.

        This method will be a no-op if the prep action does not support this
        feature.
        """
        self.prep.addInputSchema(self.process.getInputSchema())

    def getOutputNames(self, config: pexConfig.Config | None = None) -> Iterable[str]:
        """Return the names of the plots produced by this analysis tool.

        If there is a `PlotAction` defined in the produce action, these names
        will either come from the `PlotAction` if it defines a
        ``getOutputNames`` method (likely if it returns a mapping of figures),
        or a default value is used and a single figure is assumed.

        Parameters
        ----------
        config : `lsst.pex.config.Config`, optional
            Configuration of the task. This is only used if the output naming
            needs to be config-aware.

        Returns
        -------
        result : `tuple` of `str`
            Names for each plot produced by this action.
        """
        match self.produce:
            case JointAction(plot=NoPlot()):
                return tuple()
            case _HasOutputNames():
                outNames = tuple(self.produce.getOutputNames(config=config))
            case _:
                raise ValueError(f"Unsupported Action type {type(self.produce)} for getting output names")

        results = []
        suffix = self._getPlotType()
        if self.parameterizedBand:
            prefix = "_".join(x for x in ("{band}", self.identity) if x)
        else:
            prefix = f"{self.identity}" if self.identity else ""

        if outNames:
            for name in outNames:
                results.append("_".join(x for x in (prefix, name, suffix) if x))
        else:
            results.append("_".join(x for x in (prefix, suffix) if x))
        return results

    @classmethod
    def fromPipeline(
        cls,
        pipeline: str | Pipeline,
        name: str,
        fullpath: bool = False,
        instrument: Instrument | str | None = None,
    ) -> AnalysisTool | None:
        """Construct an `AnalysisTool` from a definition written in a
        `~lsst.pipe.base.Pipeline`.

        Parameters
        ----------
        pipeline : `str` or `~lsst.pipe.base.Pipeline`
            The pipeline to load the `AnalysisTool` from.
        name : `str`
            The name of the analysis tool to run. This can either be just the
            name assigned to the tool, or an absolute name in a config
            hierarchy.
        fullpath : `bool`
            Determines if the name is interpreted as an absolute path in a
            config hierarchy, or is relative to an `AnalysisTool` ``atools``
            `~lsst.pex.config.configurableActions.ConfigurableActionStructField`
            .
        instrument : `~lsst.daf.butler.instrument.Instrument` or `str` or\
            `None`
            Either a derived class object of a `lsst.daf.butler.instrument` or
            a string corresponding to a fully qualified
            `lsst.daf.butler.instrument` name or None if no instrument needs
            specified or if the pipeline contains the instrument in it.
            Defaults to None.

        Returns
        -------
        tool : `AnalysisTool`
            The loaded `AnalysisTool` as configured in the pipeline.

        Raises
        ------
        ValueError
            Raised if the config field specified does not point to an
            `AnalysisTool`.
            Raised if an instrument is specified and it conflicts with the
            pipelines instrument.
        """
        if not isinstance(pipeline, Pipeline):
            pipeline = Pipeline.fromFile(pipeline)
        # If the caller specified an instrument, verify it does not conflict
        # with the pipelines instrument, and add it to the pipeline
        if instrument is not None:
            if (pipeInstrument := pipeline.getInstrument()) and pipeInstrument != instrument:
                raise ValueError(
                    f"The supplied instrument {instrument} conflicts with the pipelines instrument "
                    f"{pipeInstrument}."
                )
            else:
                pipeline.addInstrument(instrument)
        try:
            pipelineGraph = pipeline.to_graph()
        except (FieldValidationError, ValueError) as err:
            raise ValueError(
                "There was an error instantiating the pipeline, do you need to specify an instrument?"
            ) from err
        if not fullpath:
            name = f"atools.{name}"
        for task in pipelineGraph.tasks.values():
            config = task.config
            try:
                attr = attrgetter(name)(config)
            except AttributeError:
                continue
            if not isinstance(attr, AnalysisTool):
                raise ValueError("The requested name did not refer to an analysisTool")
            return attr
        return None

    def finalize(self) -> None:
        """Run any finalization code that depends on configuration being
        complete.
        """
        pass

    def _baseFinalize(self) -> None:
        self.populatePrepFromProcess()

    def freeze(self):
        if not self.__dict__.get("_finalizeRun"):
            self.finalize()
            self.__dict__["_finalizeRun"] = True
        super().freeze()


# explicitly wrap the finalize of the base class
AnalysisTool.finalize = _finalizeWrapper(AnalysisTool.finalize, AnalysisTool)
