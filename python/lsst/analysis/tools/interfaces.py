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

from abc import abstractmethod, ABCMeta
from numbers import Number
from typing import Iterable, MutableMapping, Tuple, Type, Mapping, cast, Any
import warnings
from collections import abc

import numpy as np
from lsst.pipe.tasks.configurableActions import ConfigurableAction, ConfigurableActionField
from lsst.verify import Measurement
from matplotlib.figure import Figure
from numpy.typing import NDArray


class ScalarMeta(ABCMeta):
    def __instancecheck__(cls: ABCMeta, instance: Any) -> Any:
        return isinstance(instance, tuple(cls.mro()[1:]))


class Scalar(Number, np.number, metaclass=ScalarMeta):
    def __init__(self) -> None:
        raise NotImplementedError("Scalar is only an interface and should not be directly instantiated")


Vector = NDArray
KeyedData = MutableMapping[str, Vector | Scalar]

KeyedDataSchema = Iterable[Tuple[str, Type[Vector] | Type[Number] | type[np.number]]]


class AnalysisAction(ConfigurableAction):
    def __init_subclass__(cls, **kwargs):
        if "getInputSchema" not in dir(cls):
            raise NotImplementedError(f"Class {cls} must implement method getInputSchema")

    @abstractmethod
    def getInputSchema(self) -> KeyedDataSchema:
        raise NotImplementedError("This is not implemented on the base class")

    def getFormattedInputSchema(self, **kwargs):
        for key, typ in self.getInputSchema():
            yield key.format_map(kwargs), typ

    def addInputSchema(self, inputSchema: KeyedDataSchema) -> None:
        warnings.warn(
            f"{type(self)} does not implement adding input schemas, call will do nothing, "
            "this may be expected",
            RuntimeWarning,
        )


class KeyedDataAction(AnalysisAction):
    def __init_subclass__(cls, **kwargs):
        if "getOutputSchema" not in dir(cls):
            raise NotImplementedError(f"Class {cls} must implement method getOutputSchema")

    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        raise NotImplementedError("This is not implemented on the base class")

    @abstractmethod
    def getOutputSchema(self) -> KeyedDataSchema:
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
    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Measurement] | Measurement:
        raise NotImplementedError("This is not implemented on the base class")


class PlotAction(AnalysisAction):
    @abstractmethod
    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        raise NotImplementedError("This is not implemented on the base class")


class AnalysisTool(AnalysisAction):
    prep = ConfigurableActionField(doc="Action to run to prepare inputs", dtype=KeyedDataAction)
    process = ConfigurableActionField(
        doc="Action to process data into intended form",
    )
    post_process = ConfigurableActionField(doc="Action to perform any finalization steps")

    def __call__(
        self, data: KeyedData, **kwargs
    ) -> Mapping[str, Figure] | Figure | Mapping[str, Measurement] | Measurement:
        if (bands := kwargs.pop("bands", None)) is None:
            return self._call_single(data, **kwargs)
        results: dict[str, Any] = {}
        if self.identity is not None:
            value_key = f"{{band}}_{self.identity}"
        else:
            value_key = "{band}"
        for band in bands:
            kwargs["band"] = band
            match self._call_single(data, **kwargs):
                case abc.Mapping() as mapping:
                    results.update(mapping.items())
                case value:
                    results[value_key.format(band=band)] = value
        return results

    def _call_single(
        self, data: KeyedData, **kwargs
    ) -> Mapping[str, Figure] | Figure | Mapping[str, Measurement] | Measurement:
        self.populatePrepFromProcess()
        prepped: KeyedData = self.prep(data, **kwargs)  # type: ignore
        processed: KeyedData = self.process(prepped, **kwargs)  # type: ignore
        finalized: Mapping[str, Figure] | Figure | Mapping[str, Measurement] | Measurement =\
            self.post_process(processed, **kwargs)  # type: ignore
        return finalized

    def setDefaults(self):
        super().setDefaults()
        # imported here to avoid circular imports
        from .analysisParts.base import BasePrep, BaseProcess

        self.prep = BasePrep()
        self.process = BaseProcess()

    def getInputSchema(self) -> KeyedDataSchema:
        self.populatePrepFromProcess()
        return self.prep.getInputSchema()  # type: ignore

    def populatePrepFromProcess(self):
        cast(AnalysisAction, self.prep).addInputSchema(cast(AnalysisAction, self.process).getInputSchema())

    @classmethod
    @abstractmethod
    def getOutputDSNames(cls, **kwargs) -> Iterable[str]:
        raise NotImplementedError("This is not implemented on the base class")


class AnalysisMetric(AnalysisTool):
    post_process = ConfigurableActionField(doc="Action which returns a calculated Metric", dtype=MetricAction)

    def setDefaults(self):
        super().setDefaults()
        # imported here to avoid circular imports
        from .analysisParts.base import BaseMetricAction

        self.post_process = BaseMetricAction


class AnalysisPlot(AnalysisTool):
    post_process = ConfigurableActionField(doc="Action which returns a plot", dtype=PlotAction)
