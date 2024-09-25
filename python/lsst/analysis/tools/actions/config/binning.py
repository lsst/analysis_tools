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

__all__ = ("MagnitudeBinConfig",)

import lsst.pex.config as pexConfig


class MagnitudeBinConfig(pexConfig.Config):
    """Configuration for magnitude-binned metrics.

    Bins may only be specified in mmag widths
    """

    mag_low_min = pexConfig.Field[int](
        doc="Lower bound for the smallest magnitude bin in millimags",
        default=16000,
    )
    mag_low_max = pexConfig.Field[int](
        doc="Lower bound for the first excluded (largest) magnitude bin in millimags",
        default=31000,
    )
    mag_interval = pexConfig.Field[int](
        doc="Spacing interval of magnitude bins in millimags",
        default=1000,
        check=lambda x: x >= 10,
    )
    mag_width = pexConfig.Field[int](
        doc="Width of magnitude bins in millimags",
        default=1000,
        check=lambda x: x >= 10,
    )

    def get_bins(self):
        """Get the lower limit for each bin.

        Returns
        -------
        bins
            A list of the lower limits for each bin.
        """
        bins = list(range(self.mag_low_min, self.mag_low_max, self.mag_interval))
        return bins

    def get_name_bin(self, value_bin: int) -> str:
        """Get the default name of a bin for metrics.

        Parameters
        ----------
        value_bin
            The value of the bin (usually the lower limit).

        Returns
        -------
        name
            A string name for the bin.
        """
        prefix = str(value_bin // 1000)
        suffix = str(value_bin % 1000)
        name = f"{prefix}{'p' if suffix else ''}{suffix}"
        return name

    def validate(self):
        if not self.mag_low_max > self.mag_low_min:
            raise ValueError(f"{self.mag_low_max} !> {self.mag_low_min}")
