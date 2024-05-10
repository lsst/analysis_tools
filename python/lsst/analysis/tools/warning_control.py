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

__all__ = (
    "filterwarnings_action",
    "numpy_all_nan",
    "numpy_all_nan_slice",
    "numpy_divide_zero_divide",
    "numpy_divide_zero_log",
    "numpy_divide_zero_log10",
    "numpy_dof_zero",
    "numpy_invalid_value_cos",
    "numpy_invalid_value_divide",
    "numpy_invalid_value_log",
    "numpy_invalid_value_log10",
    "numpy_invalid_value_subtract",
    "numpy_invalid_value_sin",
    "numpy_invalid_value_sqrt",
    "numpy_mean_empty",
)

# This will be applied to all actions.
# "ignore" and "error" should be used for running and debugging, respectively.
filterwarnings_action = "ignore"

numpy_all_nan = r"All-NaN axis encountered"
numpy_all_nan_slice = r"All-NaN slice encountered"
numpy_divide_zero_divide = r"divide by zero encountered in divide"
numpy_divide_zero_log = r"divide by zero encountered in log"
numpy_divide_zero_log10 = r"divide by zero encountered in log10"
numpy_dof_zero = r"Degrees of freedom <= 0 for slice."
numpy_invalid_value_cos = r"invalid value encountered in cos"
numpy_invalid_value_divide = r"invalid value encountered in divide"
numpy_invalid_value_log = r"invalid value encountered in log"
numpy_invalid_value_log10 = r"invalid value encountered in log10"
numpy_invalid_value_subtract = r"invalid value encountered in subtract"
numpy_invalid_value_sin = r"invalid value encountered in sin"
numpy_invalid_value_sqrt = r"invalid value encountered in sqrt"
numpy_mean_empty = r"Mean of empty slice"
