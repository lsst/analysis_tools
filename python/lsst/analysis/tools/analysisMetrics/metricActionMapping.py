from __future__ import annotations

__all__ = ("MetricActionMapping",)

from collections import UserDict

from lsst.verify import Measurement


class MetricActionMapping(UserDict[str, list[Measurement]]):
    """A specialized dict for storing outputs from multiple `AnalysisMetric`
    actions.

    Keys correspond to the identifier of the action that produced the
    corresponding values. Each value is a list of `lsst.verift.Measurement`
    objects produced by that `AnalysisMetric` action.

    This object also supports the pydandic interface for serializing to and
    from json compatible objects.
    """
    def json(self):
        result = {}
        for key, value in self.items():
            result[key] = [meas.json for meas in value]
        return result

    @classmethod
    def parse_obj(cls, data):
        inst = cls()
        for key, value in data:
            inst[key] = [Measurement.deserialize(**element) for element in value]
        return inst
