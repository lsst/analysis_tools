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
r"""The context module provides the tooling to define various execution
contexts.

Execution contexts are defined by creating a subclass of `Context`. By
Convention the name of the subclass should be <Name>Context.

Multiple contexts can be joined together to create an extended Context by
or-ing them together i.e.

compoundContext = VisitContext | PlotContext

A `Context` is applied to an `AnalysisAction` by calling ``applyContext`` and
passing in a context. The `AnalysisAction` will then see if it has
corresponding method (the name of the context with the first letter lower case
i.e. visitContext). If a method is present, it will be called to activate the
context for that `AnalysisAction`. The action will then recursively walk down
through any  subsequent actions calling ``applyContext``.

When the context being applied is some compound context, a call to
``applyContext`` will happen for each `Context` individually. The variable name
of the compound `Context` is not important, and in some cases is unneeded i.e.

analysis.applyContext(VisitContext | PlotContext)

This is equivalent to:

analysis.applyContext(VisitContext)
analysis.applyContext(PlotContext)

But the ability to store multiple `Context`\ s in a variable can be useful when
linking multiple `AnalysisTool`\ s together, or creating a context that might
be frequently reused.

An `AnalysisAction` is entirely in control to what it does when it enters a
`Context` or if it even responds at all based on if it implements the context
method, and what it does inside that method.
"""

from ._baseContext import *
from ._contexts import *
