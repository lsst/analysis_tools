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

__all__ = ("StellarLocusFitAction",)

from typing import cast

import numpy as np
import smatch
from lsst.pex.config import Field
from scipy.stats import gaussian_kde

from ...interfaces import KeyedData, KeyedDataAction, KeyedDataSchema, Scalar, Vector


class depthCalculation(KeyedDataAction):
    "Determine limiting magnitude given a mag and magerr vectors"
    magKey = Field[str](doc="Key to extract magnitude vector", default="mags")
    magErrKey = Field[str](doc="Key to extract magnitude error vector", default="magErrs")
    raKey = Field[str](doc="key to extract ra vector", default="coord_ra")
    decKey = Field[str](doc="key to extract dec vector", default="coord_dec")

    signalToNoise = Field[float](
        doc="signal to nose value will compute mag where snr equals this value", default=5.0
    )
    magMin = Field[float](doc="minimum magnitude to consider for calculation", default=16.0)
    magMax = Field[float](doc="minimum magnitude to consider for calculation", default=30.0)

    def getInputSchema(self) -> KeyedDataSchema:
        return ((self.raKey, Vector), (self.decKey, Vector))

    def __call__(self, data: KeyedData, **kwargs) -> KeyedData:

        mag_snr = (2.5 / np.log(10)) / (self.signalToNoise)

        h, edges = np.histogram(data[self.magKey], bins=np.arange(self.magMin, self.magMax, 0.1))

        mag_bright_end = edges[np.argmax(h)] - 3.0

        cut = (data[self.magKey] > np.max([mag_bright_end, self.magMin])) & (data[self.magKey] < self.magMax)

        magVector = data[self.magKey][cut]
        magErrVector = data[self.magErrKey][cut]
        raVector = data[self.raKey][cut]
        decVector = data[self.decKey][cut]

        match = smatch.match_self(raVector, decVector, radius=0.2, maxmatch=1)

        delta_mag = magVector[match["i1"]] - magVector[match["i2"]]
        delta_log_magerr = np.log10(magErrVector[match["i1"]]) - np.log10(magErrVector[match["i2"]])
        # import pdb; pdb.set_trace()
        # old = np.seterr(divide='ignore',invalid='ignore')
        np.warnings.filterwarnings("ignore", r"invalid value encountered")  # type: ignore
        np.warnings.filterwarnings("ignore", r"divide by zero")  # type: ignore
        ratio = delta_log_magerr / delta_mag
        cut_nan_inf = np.isfinite(ratio) & (delta_mag > 0.5)

        # np.seterr(**old)

        if cut_nan_inf.sum() < 2:
            # logger not enough sources
            return np.nan

        kde = gaussian_kde(ratio[cut_nan_inf])

        values = np.linspace(0.0, 1.0, 1000)
        kde_values = kde.evaluate(values)
        slope = values[np.argmax(kde_values)]

        maglims = magVector - ((np.log10(magErrVector) - np.log10(mag_snr)) / slope)

        maglim = np.nanmedian(maglims)

        return cast(Scalar, maglim)
