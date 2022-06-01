from __future__ import annotations

__all__ = (
    "TabularAction",
    "ScalarAction",
    "MetricAction",
    "PlotAction",
    "NumberType",
    "Tabular",
    "Vector",
    "CalculatorAction",
    "AnalysisTool",
    "AnalysisMetric",
    "AnalysisPlot",
)

from abc import abstractmethod
from numbers import Number
from typing import Any, Iterable, Mapping, MutableMapping, NewType

import numpy as np
from lsst.pipe.tasks.configurableActions import (ConfigurableAction,
                                                 ConfigurableActionField)
from lsst.verify import Metric
from matplotlib.figure import Figure
from numpy.typing import NDArray

NumberType = Number | np.number
Vector = NewType("Vector", NDArray)
Tabular = MutableMapping[str, Vector]


class TabularAction(ConfigurableAction):
    @abstractmethod
    def getColumns(self, **kwargs) -> Iterable[str]:
        raise NotImplementedError("This is not implemented on the base class")

    @abstractmethod
    def __call__(self, table: Tabular, **kwargs) -> Tabular:
        raise NotImplementedError("This is not implemented on the base class")


class VectorAction(ConfigurableAction):
    @abstractmethod
    def getColumns(self, **kwargs) -> Iterable[str]:
        raise NotImplementedError("This is not implemented on the base class")

    @abstractmethod
    def __call__(self, table: Tabular, **kwargs) -> Vector:
        raise NotImplementedError("This is not implemented on the base class")


class ScalarAction(ConfigurableAction):
    @abstractmethod
    def getColumns(self, **kwargs) -> Iterable[str]:
        raise NotImplementedError("This is not implemented on the base class")

    @abstractmethod
    def __call__(self, table: Tabular, **kwargs) -> NumberType:
        raise NotImplementedError("This is not implemented on the base class")

    def getMask(self, **kwargs) -> Vector | slice:
        if (mask := kwargs.get("mask")) is None:
            mask = slice(None)
        return mask


CalculatorAction = TabularAction | ScalarAction | VectorAction


class MetricAction(ConfigurableAction):
    @abstractmethod
    def __call__(self, input: Tabular | NumberType, **kwargs) -> Mapping[str, Metric] | Metric:
        raise NotImplementedError("This is not implemented on the base class")


class PlotAction(ConfigurableAction):
    @abstractmethod
    def __call__(self, input: Tabular | NumberType, **kwargs) -> Mapping[str, Figure] | Figure:
        raise NotImplementedError("This is not implemented on the base class")


class AnalysisTool(ConfigurableAction):
    prep = ConfigurableActionField(doc="Action to run to prepare inputs", dtype=TabularAction)
    process = ConfigurableActionField(doc="Action to process data into intended form", dtype=CalculatorAction)
    post_process = ConfigurableActionField(doc="Action to perform any finalization steps")

    def __call__(self, table: Mapping, **kwargs) -> Any:
        prepped = self.prep(table, **kwargs)  # type: ignore
        processed = self.process(prepped, **kwargs)  # type: ignore
        finalized = self.post_process(processed, **kwargs)  # type: ignore
        return finalized

    @classmethod
    @abstractmethod
    def getOutputDSNames(cls, **kwargs) -> Iterable[str]:
        raise NotImplementedError("This is not implemented on the base class")


class AnalysisMetric(AnalysisTool):
    post_process = ConfigurableActionField(doc="Action which returns a calculated Metric", dtype=MetricAction)


class AnalysisPlot(AnalysisTool):
    post_process = ConfigurableActionField(doc="Action which returns a plot", dtype=PlotAction)


# For Demo Example #


# class AnalysisBaseConfig(Config):
#     analysisActions = ConfigurableActionStructField(doc="AnalysisActions to run with this Task")
#
#     def setDefaults(self):
#         self.analysisActions.te1 = TE1
#
#
# class AnalysisPipelineTask(PipelineTask):
#     def run(self, table, **kwargs) -> Struct:
#         results = Struct()
#         for name, action in self.config.analysisActions:  # type: ignore
#             match action(table, **kwargs):
#                 case Mapping(val):
#                     for n, v in val.items():
#                         setattr(results, n, v)
#                 case value:
#                     setattr(results, name, value)
#         return results
