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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from __future__ import annotations

import warnings
from typing import cast
from unittest import TestCase, main

import astropy.units as apu
import numpy as np

import lsst.utils.tests
from lsst.analysis.tools.actions.scalar import MeanAction, MedianAction
from lsst.analysis.tools.contexts import Context
from lsst.analysis.tools.interfaces import (
    AnalysisAction,
    AnalysisTool,
    KeyedData,
    KeyedDataAction,
    KeyedDataSchema,
    KeyedResults,
    Scalar,
    ScalarAction,
    Vector,
)
from lsst.pex.config import Field
from lsst.pex.config.configurableActions import ConfigurableActionField, ConfigurableActionStructField
from lsst.verify import Measurement


class MedianContext(Context):
    """Test Context for median"""

    pass


class MeanContext(Context):
    """Test Context for mean"""

    pass


class MultiplyContext(Context):
    """Test Context to multiply results"""

    pass


class DivideContext(Context):
    """Test Context to divide result"""

    pass


class TestAction1(KeyedDataAction):
    __test__ = False  # Tell pytest that this is *not* a test suite

    multiple = ConfigurableActionStructField[ScalarAction](doc="Multiple Actions")

    def getInputSchema(self) -> KeyedDataSchema:
        return (("a", Vector),)

    def medianContext(self) -> None:
        self.multiple.a = MedianAction(vectorKey="a")
        self.multiple.b = MedianAction(vectorKey="a")
        self.multiple.c = MedianAction(vectorKey="a")

    def meanContext(self) -> None:
        self.multiple.a = MeanAction(vectorKey="a")
        self.multiple.b = MeanAction(vectorKey="a")

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result = np.array([action(data, **kwargs) for action in self.multiple])
        return cast(KeyedData, {"b": result})


class TestAction2(KeyedDataAction):
    __test__ = False  # Tell pytest that this is *not* a test suite

    single = ConfigurableActionField[ScalarAction](doc="Single Action")
    multiplier = ConfigurableActionField[AnalysisAction](doc="Multiplier")

    def getInputSchema(self) -> KeyedDataSchema:
        return (("b", Vector),)

    def setDefaults(self) -> None:
        super().setDefaults()
        # This action remains constant in every context, but it will be
        # recursively configured with any contexts
        self.multiplier = TestAction3()

    def medianContext(self) -> None:
        self.single = MedianAction(vectorKey="b")

    def meanContext(self) -> None:
        self.single = MeanAction(vectorKey="b")

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:
        result = self.single(data, **kwargs)
        intermediate = self.multiplier(cast(KeyedData, {"c": result}), **kwargs)["d"]
        intermediate = cast(Measurement, intermediate)
        assert intermediate.quantity is not None
        result *= cast(Scalar, intermediate.quantity.value)
        return {"c": result}


class TestAction3(KeyedDataAction):
    __test__ = False  # Tell pytest that this is *not* a test suite

    multiplier = Field[float](doc="Number to multiply result by")

    def getInputSchema(self) -> KeyedDataSchema:
        return (("c", Vector),)

    def multiplyContext(self):
        self.multiplier = 2

    def divideContext(self):
        self.multiplier = 0.5

    def __call__(self, data: KeyedData, **kwargs) -> KeyedResults:
        result = Measurement("TestMeasurement", cast(Scalar, data["c"]) * self.multiplier * apu.Unit("count"))
        return {"d": result}


class TestAnalysisTool(AnalysisTool):
    __test__ = False  # Tell pytest that this is *not* a test suite

    def setDefaults(self) -> None:
        self.prep = TestAction1()
        self.process = TestAction2()
        self.produce = TestAction3()


class ContextTestCase(TestCase):
    def setUp(self) -> None:
        super().setUp()
        self.array = np.arange(20)
        self.input = cast(KeyedData, {"a": self.array})

    def testContext1(self):
        tester = TestAnalysisTool()
        # test applying contexts serially
        tester.applyContext(MultiplyContext())

        # verify assignment syntax works to support Yaml
        # normally this should be called as a function in python
        tester.applyContext = MedianContext
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = cast(Measurement, tester(self.input)["d"])
            assert result.quantity is not None
        self.assertEqual(result.quantity.value, 361)

    def testContext2(self):
        tester2 = TestAnalysisTool()
        compound = MeanContext | DivideContext
        tester2.applyContext(compound)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = cast(Measurement, tester2(self.input)["d"])
            assert result.quantity is not None
        self.assertEqual(result.quantity.value, 22.5625)


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    main()
