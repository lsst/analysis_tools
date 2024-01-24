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

__all__ = ("MatrixPlot",)

from typing import TYPE_CHECKING, Any, Mapping

import astropy.visualization as apViz
import matplotlib.patheffects as mpl_path_effects
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization.mpl_normalize import ImageNormalize
from lsst.pex.config import ChoiceField, Config, ConfigDictField, DictField, Field, ListField

from ...interfaces import PlotAction, Vector
from .plotUtils import addPlotInfo

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from ...interfaces import KeyedData, KeyedDataSchema


class GuideLinesConfig(Config):
    lines = DictField[float, str](
        doc=("Dictionary of x/y-values and the labels where vertical/horizontal lines are drawn."),
        optional=False,
    )

    color = Field[str](
        doc="The color of the lines and labels.",
        default="red",
    )

    outlineColor = Field[str](
        doc="The color of the outline around the lines and labels.",
        default="white",
    )

    linestyle = Field[str](
        doc="The style of the lines.",
        default="--",
    )


class MatrixPlot(PlotAction):
    """Make the plot of a matrix (2D array).

    Notes
    -----
    The `xAxisTickLabels` and `yAxisTickLabels` attributes of this class serve
    as dictionaries to map axis tick positions to their corresponding labels.
    If any positions do not align with major ticks (either provided by
    `x/yAxisTickValues` or automatically set by matplotlib), they will be
    designated as minor ticks. Thus, these tick labels operate independently,
    meaning their corresponding positions do not need to match those in
    `x/yAxisTickValues` or anything else. The code automatically adjusts to
    handle any overlaps caused by user input and across various plotting
    scenarios.
    Note that when `component1Key` and `component2Key` are specified, the x and
    y tick values and labels will be dynamically configured, thereby
    eliminating the need for providing `x/yAxisTickValues` and
    `x/yAxisTickLabels`.
    """

    inputDim = ChoiceField[int](
        doc="The dimensionality of the input data.",
        default=1,
        allowed={
            1: "1D inputs are automatically reshaped into square 2D matrices.",
            2: "2D inputs are directly utilized as is.",
        },
        optional=True,
    )

    matrixKey = Field[str](
        doc="The key for the input matrix.",
        default="matrix",
    )

    matrixOrigin = ChoiceField[str](
        doc="Determines the starting corner ('upper', 'lower') for matrix plots.",
        default="upper",
        allowed={
            "upper": "The origin is at the upper left corner.",
            "lower": "The origin is at the lower left corner.",
        },
        optional=True,
    )

    component1Key = Field[str](
        doc="The key to access a list of names for the first component set in a correlation analysis. This "
        "will be used to determine x-axis tick values and labels.",
        default=None,
        optional=True,
    )

    component2Key = Field[str](
        doc="The key to access a list of names for the second component set in a correlation analysis. This "
        "will be used to determine y-axis tick values and labels.",
    )

    xAxisLabel = Field[str](
        doc="The label to use for the x-axis.",
        default="",
        optional=True,
    )

    yAxisLabel = Field[str](
        doc="The label to use for the y-axis.",
        default="",
        optional=True,
    )

    axisLabelFontSize = Field[float](
        doc="The font size for the axis labels.",
        default=9,
        optional=True,
    )

    colorbarLabel = Field[str](
        doc="The label to use for the colorbar.",
        default="",
        optional=True,
    )

    colorbarLabelFontSize = Field[float](
        doc="The font size for the colorbar label.",
        default=10,
        optional=True,
    )

    colorbarTickLabelFontSize = Field[float](
        doc="The font size for the colorbar tick labels.",
        default=8,
        optional=True,
    )

    colorbarCmap = ChoiceField[str](
        doc="The colormap to use for the colorbar.",
        default="viridis",
        allowed={name: name for name in plt.colormaps()},
        optional=True,
    )

    vmin = Field[float](
        doc="The vmin value for the colorbar.",
        default=None,
        optional=True,
    )

    vmax = Field[float](
        doc="The vmax value for the colorbar.",
        default=None,
        optional=True,
    )

    figsize = ListField[float](
        doc="The size of the figure.",
        default=[5, 5],
        maxLength=2,
        optional=True,
    )

    title = Field[str](
        doc="The title of the figure.",
        default="",
        optional=True,
    )

    titleFontSize = Field[float](
        doc="The font size for the title.",
        default=10,
        optional=True,
    )

    xAxisTickValues = ListField[float](
        doc="List of x-axis tick values. If not set, the ticks will be set automatically by matplotlib.",
        default=None,
        optional=True,
    )

    xAxisTickLabels = DictField[float, str](
        doc="Dictionary mapping x-axis tick positions to their corresponding labels. For behavior details, "
        "refer to the 'Notes' section of the class docstring.",
        default=None,
        optional=True,
    )

    yAxisTickValues = ListField[float](
        doc="List of y-axis tick values. If not set, the ticks will be set automatically by matplotlib.",
        default=None,
        optional=True,
    )

    yAxisTickLabels = DictField[float, str](
        doc="Dictionary mapping y-axis tick positions to their corresponding labels. For behavior details, "
        "refer to the 'Notes' section of the class docstring.",
        default=None,
        optional=True,
    )

    tickLabelsFontSize = Field[float](
        doc="The font size for the tick labels.",
        default=8,
        optional=True,
    )

    tickLabelsRotation = Field[float](
        doc="The rotation of the tick labels.",
        default=0,
        optional=True,
    )

    setPositionsAtPixelBoundaries = Field[bool](
        doc="Whether to consider the positions at the pixel boundaries rather than the center of the pixel.",
        default=False,
        optional=True,
    )

    hideMajorTicks = ListField[str](
        doc="List of axis names for which to hide the major ticks. The options to include in the list are "
        "'x' and 'y'. This does not affect the visibility of major tick 'labels'. For example, setting this "
        "field to ['x', 'y'] will hide both major ticks.",
        default=[],
        maxLength=2,
        itemCheck=lambda s: s in ["x", "y"],
        optional=True,
    )

    hideMinorTicks = ListField[str](
        doc="List of axis names for which to hide the minor ticks. The options to include in the list are "
        "'x' and 'y'. This does not affect the visibility of minor tick labels. For example, setting this "
        "field to ['x', 'y'] will hide both minor ticks.",
        default=[],
        maxLength=2,
        itemCheck=lambda s: s in ["x", "y"],
        optional=True,
    )

    dpi = Field[int](
        doc="The resolution of the figure.",
        default=300,
        optional=True,
    )

    guideLines = ConfigDictField[str, GuideLinesConfig](
        doc="Dictionary of guide lines for the x and y axes. The keys are 'x' and 'y', and the values are "
        "instances of `GuideLinesConfig`.",
        default={},
        dictCheck=lambda d: all([k in ["x", "y"] for k in d]),
        optional=True,
    )

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector]]] = []
        base.append((self.matrixKey, Vector))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData, **kwargs: Any) -> None:
        # Check that the input data contains all the required keys.
        needed = set(k[0] for k in self.getInputSchema())
        if not needed.issubset(data.keys()):
            raise ValueError(f"Input data does not contain all required keys: {self.getInputSchema()}")
        # Check the input data is a matrix, i.e. a 2d array.
        if not isinstance(data[self.matrixKey], np.ndarray) and data[self.matrixKey].ndim != 2:
            raise ValueError(f"Input data is not a 2d array: {data[self.matrixKey]}")
        # Check that the keyword arguments are valid.
        acceptableKwargs = {"plotInfo", "skymap", "band", "metric_tags", "fig"}
        if not set(kwargs).issubset(acceptableKwargs):
            raise ValueError(
                f"Only the following keyword arguments are allowed: {acceptableKwargs}. Got: {kwargs}"
            )
        # Check that if one component key is provided, the other must be too.
        if (self.component1Key is not None and self.component2Key is None) or (
            self.component1Key is None and self.component2Key is not None
        ):
            raise ValueError(
                "Both 'component1Key' and 'component2Key' must be provided together if either is provided."
            )
        # Check that if component keys are provided, any of the tick values or
        # labels are not and vice versa.
        if (self.component1Key is not None and self.component2Key is not None) and (
            self.xAxisTickValues is not None
            or self.yAxisTickValues is not None
            or self.xAxisTickLabels is not None
            or self.yAxisTickLabels is not None
        ):
            raise ValueError(
                "If 'component1Key' and 'component2Key' are provided, 'xAxisTickValues', "
                "'yAxisTickValues', 'xAxisTickLabels', and 'yAxisTickLabels' should not be "
                "provided as they will be dynamically configured."
            )

    def makePlot(self, data: KeyedData, plotInfo: Mapping[str, str] | None = None, **kwargs: Any) -> Figure:
        """
        Plot a matrix of values.

        Parameters
        ----------
        data : `~lsst.analysis.tools.interfaces.KeyedData`
            The data to plot.
        plotInfo : `dict`, optional
            A dictionary of information about the data being plotted.
        **kwargs
            Additional keyword arguments to pass to the plot.

        Returns
        -------
        fig : `~matplotlib.figure.Figure`
            The resulting figure.
        """
        # Retrieve the matrix info from the input data.
        matrix = data[self.matrixKey]

        # Fetch the components between which the correlation is calculated.
        if self.component1Key is not None and self.component2Key is not None:
            comp1 = data[self.component1Key]
            comp2 = data[self.component2Key]

        if self.inputDim == 1:
            # Calculate the size of the square.
            square_size = int(np.sqrt(matrix.size))
            # Reshape into a square array.
            matrix = matrix.reshape(square_size, square_size)
            if self.component1Key is not None and self.component2Key is not None:
                comp1 = comp1.reshape(square_size, square_size)
                comp2 = comp2.reshape(square_size, square_size)

        # Calculate default limits only if needed.
        if self.vmin is None or self.vmax is None:
            default_limits = apViz.PercentileInterval(98.0).get_limits(np.abs(matrix.flatten()))
        else:
            default_limits = (None, None)

        # Set the value range using overrides or defaults.
        vrange = (
            default_limits[0] if self.vmin is None else self.vmin,
            default_limits[1] if self.vmax is None else self.vmax,
        )

        # Allow for the figure object to be passed in.
        fig = kwargs.get("fig")
        if fig is None:
            fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
            ax = fig.add_subplot(111)
        else:
            ax = fig.gca()

        if self.title:
            ax.set_title(self.title, fontsize=self.titleFontSize)

        if self.xAxisLabel:
            ax.set_xlabel(self.xAxisLabel, fontsize=self.axisLabelFontSize)

        if self.yAxisLabel:
            ax.set_ylabel(self.yAxisLabel, fontsize=self.axisLabelFontSize)

        # Set the colorbar and draw the image.
        norm = ImageNormalize(vmin=vrange[0], vmax=vrange[1])
        img = ax.imshow(
            matrix, interpolation="none", norm=norm, origin=self.matrixOrigin, cmap=self.colorbarCmap
        )

        # Calculate the aspect ratio of the image.
        ratio = matrix.shape[0] / matrix.shape[1]

        # Add the colorbar flush with the image axis.
        cbar = fig.colorbar(img, fraction=0.0457 * ratio, pad=0.04)

        # Set the colorbar label and its font size.
        cbar.set_label(self.colorbarLabel, fontsize=self.colorbarLabelFontSize)

        # Set the colorbar tick label font size.
        cbar.ax.tick_params(labelsize=self.colorbarTickLabelFontSize)

        # If requested, we shift all the positions by 0.5 considering the
        # zero-point at a pixel boundary rather than the center of the pixel.
        shift = 0.5 if self.setPositionsAtPixelBoundaries else 0

        if self.component1Key is not None and self.component2Key is not None:
            xAxisTickValues = np.arange(matrix.shape[0] + shift)
            yAxisTickValues = np.arange(matrix.shape[1] + shift)
            xAxisTickLabels = {key + shift: str(val) for key, val in zip(range(matrix.shape[0]), comp1[0, :])}
            yAxisTickLabels = {key + shift: str(val) for key, val in zip(range(matrix.shape[1]), comp2[:, 0])}
        else:
            xAxisTickValues = self.xAxisTickValues
            yAxisTickValues = self.yAxisTickValues
            xAxisTickLabels = self.xAxisTickLabels
            yAxisTickLabels = self.yAxisTickLabels

        # If the tick values are not provided, retrieve them from the axes.
        xticks = xAxisTickValues if xAxisTickValues is not None else ax.xaxis.get_ticklocs()
        yticks = yAxisTickValues if yAxisTickValues is not None else ax.yaxis.get_ticklocs()

        # Retrieve the current limits of the x and y axes.
        xlim, ylim = ax.get_xlim(), ax.get_ylim()

        # Filter out tick locations that fall outside the current x/y-axis
        # limits to ensures that only tick locations within the visible range
        # are kept.
        xticks = np.array([tick for tick in xticks if min(xlim) <= tick - shift <= max(xlim)])
        yticks = np.array([tick for tick in yticks if min(ylim) <= tick - shift <= max(ylim)])
        tick_data = {
            "x": (
                xticks - shift,
                np.array(list(xAxisTickLabels.keys())) - shift if xAxisTickLabels else None,
                list(xAxisTickLabels.values()) if xAxisTickLabels else None,
            ),
            "y": (
                yticks - shift,
                np.array(list(yAxisTickLabels.keys())) - shift if yAxisTickLabels else None,
                list(yAxisTickLabels.values()) if yAxisTickLabels else None,
            ),
        }

        for dim, axis in [("x", ax.xaxis), ("y", ax.yaxis)]:
            # Get the major tick positions and labels.
            major_tick_values, positions, labels = tick_data[dim]

            # Set major ticks.
            axis.set_ticks(major_tick_values, minor=False)

            # Set tick labels while compensating for the potential shift in the
            # tick positions and removing trailing zeros and the decimal point
            # for integer values.
            axis.set_ticklabels(
                [
                    f"{tick + shift:.0f}" if (tick + shift).is_integer() else f"{tick + shift}"
                    for tick in axis.get_ticklocs()
                ],
                fontsize=self.tickLabelsFontSize,
            )

            # Check if positions are provided.
            if positions is not None:
                # Assign specified positions as minor ticks.
                axis.set_ticks(positions, minor=True)

                # Conditionally assign labels to major and/or minor ticks.
                if labels is not None:
                    # Create a lookup for positions to labels.
                    positions_labels_lookup = {
                        p: l if p in major_tick_values else "" for p, l in zip(positions, labels)
                    }
                    # Generate labels for major ticks, leaving blanks for
                    # non-major positions.
                    major_labels = [
                        "" if m not in positions_labels_lookup else positions_labels_lookup[m]
                        for m in major_tick_values
                    ]
                    # Generate labels for minor ticks, excluding those
                    # designated as major.
                    minor_labels = ["" if p in major_tick_values else l for p, l in zip(positions, labels)]

                    # Apply labels to major ticks if any exist.
                    if any(e for e in major_labels if e):
                        axis.set_ticklabels(major_labels, minor=False, fontsize=self.tickLabelsFontSize)
                    else:
                        # If no major labels, clear major tick labels.
                        axis.set_ticklabels("")

                    # Apply labels to minor ticks if any exist.
                    if any(e for e in minor_labels if e):
                        axis.set_ticklabels(minor_labels, minor=True, fontsize=self.tickLabelsFontSize)

            if dim in self.hideMajorTicks:
                # Remove major tick marks for asthetic reasons.
                axis.set_tick_params(which="major", length=0)

            if dim in self.hideMinorTicks:
                # Remove minor tick marks for asthetic reasons.
                axis.set_tick_params(which="minor", length=0)

            # Rotate the tick labels by the specified angle.
            ax.tick_params(axis=dim, rotation=self.tickLabelsRotation)

        # Add vertical and horizontal lines if provided.
        if "x" in self.guideLines:
            xLines = self.guideLines["x"]
            for x, label in xLines.lines.items():
                ax.axvline(x=x - shift, color=xLines.outlineColor, linewidth=2, alpha=0.6)
                ax.axvline(
                    x=x - shift, color=xLines.color, linestyle=xLines.linestyle, linewidth=1, alpha=0.85
                )
                label = ax.text(
                    x - shift,
                    0.03,
                    label,
                    rotation=90,
                    color=xLines.color,
                    transform=ax.get_xaxis_transform(),
                    horizontalalignment="right",
                    alpha=0.9,
                )
                # Add a distinct outline around the label for better visibility
                # in various backgrounds.
                label.set_path_effects(
                    [
                        mpl_path_effects.Stroke(linewidth=2, foreground=xLines.outlineColor, alpha=0.8),
                        mpl_path_effects.Normal(),
                    ]
                )

        if "y" in self.guideLines:
            yLines = self.guideLines["y"]
            for y, label in yLines.lines.items():
                ax.axhline(y=y - shift, color=yLines.outlineColor, linewidth=2, alpha=0.6)
                ax.axhline(
                    y=y - shift, color=yLines.color, linestyle=yLines.linestyle, linewidth=1, alpha=0.85
                )
                label = ax.text(
                    0.03,
                    y - shift,
                    label,
                    color=yLines.color,
                    transform=ax.get_yaxis_transform(),
                    verticalalignment="bottom",
                    alpha=0.9,
                )
                # Add a distinct outline around the label for better visibility
                # in various backgrounds.
                label.set_path_effects(
                    [
                        mpl_path_effects.Stroke(linewidth=2, foreground=yLines.outlineColor, alpha=0.8),
                        mpl_path_effects.Normal(),
                    ]
                )
        # Add plot info if provided.
        if plotInfo is not None:
            fig = addPlotInfo(fig, plotInfo)

        return fig
