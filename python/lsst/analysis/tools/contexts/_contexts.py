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

"""This is a module where concrete Contexts should be defined. These should
be a subclass of `Context`, and should contain a description of what the
context is for as it's docstring.
"""
__all__ = (
    "VisitContext",
    "CoaddContext",
)

from ._baseContext import Context


class VisitContext(Context):
    """A context which indicates `AnalysisAction`s are being run in the context
    of visit level data.
    """

    pass


class CoaddContext(Context):
    """A context which indicates `AnalysisAction`s are being run in the context
    of coadd level data.
    """

    pass
