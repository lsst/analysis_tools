from __future__ import annotations

__all__ = (
    "AnalysisAction",
    "KeyedDataAction",
    "ScalarAction",
    "MetricAction",
    "PlotAction",
    "Scalar",
    "KeyedData",
    "KeyedDataSchema",
    "Vector",
    "AnalysisTool",
    "AnalysisMetric",
    "AnalysisPlot",
)

from abc import abstractmethod
from numbers import Number
from typing import Iterable, MutableMapping, Tuple, Type, Mapping, Any

import numpy as np
from lsst.pipe.tasks.configurableActions import ConfigurableAction, ConfigurableActionField
from lsst.verify import Measurement
from matplotlib.figure import Figure
from numpy.typing import NDArray

Scalar = Number | np.number
Vector = NDArray
KeyedData = MutableMapping[str, Vector | Scalar]

KeyedDataSchema = Iterable[Tuple[str, Type[Vector] | Type[Scalar]]]


class AnalysisAction(ConfigurableAction):
    @abstractmethod
    def getInputSchema(self, **kwargs) -> KeyedDataSchema:
        raise NotImplementedError("This is not implemented on the base class")


class KeyedDataAction(AnalysisAction):
    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        raise NotImplementedError("This is not implemented on the base class")


class VectorAction(AnalysisAction):
    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> Vector:
        raise NotImplementedError("This is not implemented on the base class")


class ScalarAction(AnalysisAction):
    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> Scalar:
        raise NotImplementedError("This is not implemented on the base class")

    def getMask(self, **kwargs) -> Vector | slice:
        if (mask := kwargs.get("mask")) is None:
            mask = slice(None)
        return mask


class MetricAction(AnalysisAction):
    @abstractmethod
    def __call__(self, input: KeyedData | Scalar, **kwargs) -> Mapping[str, Measurement] | Measurement:
        raise NotImplementedError("This is not implemented on the base class")


class PlotAction(AnalysisAction):
    @abstractmethod
    def __call__(self, input: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        raise NotImplementedError("This is not implemented on the base class")


class AnalysisTool(ConfigurableAction):
    prep = ConfigurableActionField(doc="Action to run to prepare inputs", dtype=KeyedDataAction)
    process = ConfigurableActionField(
        doc="Action to process data into intended form",
    )
    post_process = ConfigurableActionField(doc="Action to perform any finalization steps")

    def __call__(self, data: KeyedData, **kwargs) -> Any:
        prepped = self.prep(data, **kwargs)  # type: ignore
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
