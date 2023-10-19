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

__all__ = (
    "CpAnalysisConfig",
    "CpAnalysisCalibConfig",
    "CpBiasAnalysisTask",
    "CpDarkAnalysisTask",
    "CpFlatAnalysisTask",
    "CpPtcAnalysisTask",
)

import numpy as np
from lsst.pipe.base import connectionTypes as cT

from ..interfaces import AnalysisBaseConfig, AnalysisBaseConnections, AnalysisPipelineTask


class CpAnalysisConnections(
    AnalysisBaseConnections,
    dimensions=("instrument",),
    defaultTemplates={"outputName": "cpTest1"},
):
    data = cT.Input(
        doc="Input data from cp_verify.",
        name="verifyCalibStats",
        storageClass="StructuredDataDict",
        dimensions=("instrument",),
    )
    rawData = cT.Input(
        doc="Raw values that should be collated earlier.",
        name="verifyCalibDetStats",
        storageClass="StructuredDataDict",
        dimensions=(
            "instrument",
            "exposure",
            "detector",
        ),
        deferLoad=True,
        multiple=True,
    )

    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class CpAnalysisConfig(AnalysisBaseConfig, pipelineConnections=CpAnalysisConnections):
    def setDefaults(self):
        super().setDefaults()


class CpAnalysisCalibConnections(
    AnalysisBaseConnections,
    dimensions=("instrument",),
    defaultTemplates={"outputName": "cpCalibTest"},
):
    data = cT.Input(
        doc="Input data from cp_verify.",
        name="verifyCalibStats",
        storageClass="StructuredDataDict",
        dimensions=("instrument",),
    )
    rawData = cT.Input(
        doc="Raw values that should be collated earlier.",
        name="verifyCalibDetStats",
        storageClass="StructuredDataDict",
        dimensions=(
            "instrument",
            "detector",
        ),
        deferLoad=True,
        multiple=True,
    )

    camera = cT.PrerequisiteInput(
        doc="Input camera to use for focal plane geometry.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)


class CpAnalysisCalibConfig(AnalysisBaseConfig, pipelineConnections=CpAnalysisCalibConnections):
    def setDefaults(self):
        super().setDefaults()


class CpAnalysisTask(AnalysisPipelineTask):
    ConfigClass = CpAnalysisConfig
    _DefaultName = "cpAnalysisTask"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # Docs inherited from base class.
        inputs = butlerQC.get(inputRefs)
        dataId = butlerQC.quantum.dataId

        data = self.loadData(inputs["rawData"])
        # import pdb; pdb.set_trace()

        band = None
        if "physical_filter" in dataId.keys():
            band = dataId["physical_filter"]
        instrument = dataId["instrument"]
        plotInfo = {
            "run": inputs["rawData"][0].ref.run,
            "tableName": f"{instrument}: bias",
            "bands": [band],
        }
        outputs = self.run(data=data, plotInfo=plotInfo, camera=inputs["camera"])
        butlerQC.put(outputs, outputRefs)

    def loadData(self, rawData):
        """Convert cp_verify format to analysis_tools format.

        Parameters
        ----------
        rawData : `dict` [ `dict` [ `dict` [ `dict` ]]]
            The nested dictionaries holding the individual
            measurements.

        Returns
        -------
        repackedData : `KeyedData`
            Repacked data.
        """
        raise NotImplementedError("Subclasses must implement data loading.")


class CpBiasAnalysisTask(CpAnalysisTask):
    ConfigClass = CpAnalysisConfig
    _DefaultName = "cpBiasAnalysisTask"

    def loadData(self, rawData):
        """Convert cp_verify format to analysis_tools format.

        Parameters
        ----------
        rawData : `dict` [ `dict` [ `dict` [ `dict` ]]]
            The nested dictionaries holding the individual
            measurements.

        Returns
        -------
        repackedData : `KeyedData`
            Repacked data.
        """
        # Probably should define a schema. For bias, it is:
        repackedData = {
            "detector": [],
            "amplifier": [],
            "MEAN": [],
            "NOISE": [],
            "CR_NOISE": [],
            "READ_NOISE": [],
            "VERIFY_MEAN": [],
            "VERIFY_NOISE": [],
            "VERIFY_CR_NOISE": [],
            "VERIFY_RN": [],
        }
        for entry in rawData:
            data = entry.get()
            detector = entry.dataId["detector"]

            # These keys must exist
            ampData = data["AMP"]
            mdData = data["METADATA"]
            verifyData = data["VERIFY"]

            for ampName in ampData:
                repackedData["detector"].append(detector)
                repackedData["amplifier"].append(ampName)

                for key in ("MEAN", "NOISE", "CR_NOISE"):
                    repackedData[key].append(ampData[ampName][key])
                    # These are booleans.
                    repackedData[f"VERIFY_{key}"].append(int(verifyData["AMP"][ampName][key]))
                repackedData["VERIFY_RN"].append(int(verifyData["AMP"][ampName]["READ_NOISE_CONSISTENT"]))
                repackedData["READ_NOISE"].append(mdData["RESIDUAL STDEV"][ampName])

        repackedData["statMask"] = np.ones_like(repackedData["detector"])
        return repackedData


class CpDarkAnalysisTask(CpAnalysisTask):
    ConfigClass = CpAnalysisConfig
    _DefaultName = "cpDarkAnalysisTask"

    def loadData(self, rawData):
        """Convert cp_verify format to analysis_tools format.

        Parameters
        ----------
        rawData : `dict` [ `dict` [ `dict` [ `dict` ]]]
            The nested dictionaries holding the individual
            measurements.

        Returns
        -------
        repackedData : `KeyedData`
            Repacked data.
        """
        # Probably should define a schema. For bias, it is:
        repackedData = {
            "detector": [],
            "amplifier": [],
            "MEAN": [],
            "NOISE": [],
            "CR_NOISE": [],
            "READ_NOISE": [],
            "VERIFY_MEAN": [],
            "VERIFY_NOISE": [],
            "VERIFY_CR_NOISE": [],
            "VERIFY_RN": [],
        }
        for entry in rawData:
            data = entry.get()
            detector = entry.dataId["detector"]

            # These keys must exist
            ampData = data["AMP"]
            mdData = data["METADATA"]
            verifyData = data["VERIFY"]

            for ampName in ampData:
                repackedData["detector"].append(detector)
                repackedData["amplifier"].append(ampName)

                for key in ("MEAN", "NOISE", "CR_NOISE"):
                    repackedData[key].append(ampData[ampName][key])
                    # These are booleans.
                    repackedData[f"VERIFY_{key}"].append(int(verifyData["AMP"][ampName][key]))
                repackedData["VERIFY_RN"].append(int(verifyData["AMP"][ampName]["READ_NOISE_CONSISTENT"]))
                repackedData["READ_NOISE"].append(mdData["RESIDUAL STDEV"][ampName])

        repackedData["statMask"] = np.ones_like(repackedData["detector"])
        return repackedData


class CpFlatAnalysisTask(CpAnalysisTask):
    ConfigClass = CpAnalysisConfig
    _DefaultName = "cpFlatAnalysisTask"

    def loadData(self, rawData):
        """Convert cp_verify format to analysis_tools format.

        Parameters
        ----------
        rawData : `dict` [ `dict` [ `dict` [ `dict` ]]]
            The nested dictionaries holding the individual
            measurements.

        Returns
        -------
        repackedData : `KeyedData`
            Repacked data.
        """
        # Probably should define a schema. For bias, it is:
        repackedData = {
            "detector": [],
            "amplifier": [],
            "MEAN": [],
            "NOISE": [],
            "VERIFY_NOISE": [],
        }
        for entry in rawData:
            data = entry.get()
            detector = entry.dataId["detector"]

            # These keys must exist
            ampData = data["AMP"]
            verifyData = data["VERIFY"]

            for ampName in ampData:
                repackedData["detector"].append(detector)
                repackedData["amplifier"].append(ampName)

                for key in ("NOISE",):
                    repackedData[key].append(ampData[ampName][key])
                    # These are booleans.
                    repackedData[f"VERIFY_{key}"].append(int(verifyData["AMP"][ampName][key]))

        repackedData["statMask"] = np.ones_like(repackedData["detector"])
        return repackedData


class CpPtcAnalysisTask(CpAnalysisTask):
    ConfigClass = CpAnalysisCalibConfig
    _DefaultName = "cpPtcAnalysisTask"

    def loadData(self, rawData):
        """Convert cp_verify format to analysis_tools format.

        Parameters
        ----------
        rawData : `dict` [ `dict` [ `dict` [ `dict` ]]]
            The nested dictionaries holding the individual
            measurements.

        Returns
        -------
        repackedData : `KeyedData`
            Repacked data.
        """
        # Probably should define a schema. For bias, it is:
        repackedData = {
            "detector": [],
            "amplifier": [],
            "PTC_GAIN": [],
            "PTC_BFE_A00": [],
            "PTC_NOISE": [],
            "PTC_TURNOFF": [],
        }
        for entry in rawData:
            data = entry.get()
            detector = entry.dataId["detector"]

            # These keys must exist
            ampData = data["AMP"]
            # Ignore verifyData for now.
            # verifyData = data["VERIFY"]

            for ampName in ampData:
                repackedData["detector"].append(detector)
                repackedData["amplifier"].append(ampName)

                for key in ("PTC_GAIN", "PTC_BFE_A00", "PTC_NOISE", "PTC_TURNOFF"):
                    repackedData[key].append(ampData[ampName][key])

        repackedData["statMask"] = np.ones_like(repackedData["detector"])
        return repackedData
