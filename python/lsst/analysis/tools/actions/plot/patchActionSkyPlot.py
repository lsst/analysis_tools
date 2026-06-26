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
    "PerPatchMetricConfig",
    "PerPatchPropertyMapPlot",
)

import json
import logging
import math
from collections.abc import Iterable, Mapping
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import skyproj
from matplotlib.figure import Figure

import lsst.pex.config as pexConfig
import lsst.sphgeom as sphgeom
from lsst.analysis.tools.interfaces import KeyedData, KeyedDataSchema, PlotAction, Vector
from lsst.utils.plotting import make_figure

from .plotUtils import addPlotInfo

_LOG = logging.getLogger(__name__)


class PerPatchMetricConfig(pexConfig.Config):

    description = pexConfig.Field[str](doc="Description of the metric to use as e.g. a plot title")
    key = pexConfig.Field[str](doc="Key of the metric plot")
    vmin = pexConfig.RangeField[float](
        doc="Minimum value to plot for the metric",
        min=-np.inf,
        max=np.inf,
        inclusiveMin=False,
        inclusiveMax=False,
    )
    vmax = pexConfig.RangeField[float](
        doc="Maximum value to plot for the metric",
        min=-np.inf,
        max=np.inf,
        inclusiveMin=False,
        inclusiveMax=False,
    )

    def validate(self):
        super().validate()
        if not self.vmax > self.vmin:
            raise pexConfig.FieldValidationError(f"{self.vmin=} must be greater than {self.vmax=}")


class PerPatchPropertyMapPlot(PlotAction):
    """Make per-patch sky plots for multiple metrics.

    The input data can span multiple tracts but should not cover more than
    a few thousand patches, or else the plot will be too crowded."""

    metrics = pexConfig.ConfigDictField[str, PerPatchMetricConfig](
        doc="Configuration for the metrics to plot",
    )
    cmap_metric = pexConfig.Field[str](
        doc="Name of the color map for the per-patch metric",
        default="gray",
    )
    dec_max = pexConfig.Field[float](doc="Maximum dec for the plot in degrees", optional=True)
    dec_min = pexConfig.Field[float](doc="Minimum dec for the plot in degrees", optional=True)

    dynamicOutputNames: bool = True

    figure_colorbar_width = pexConfig.Field[float](
        doc="Fractional width of the colorbar",
        default=0.06,
        check=lambda x: 0 < x < 1,
    )
    figure_dpi = pexConfig.Field[float](
        doc="Dots per inch resolution of the image",
        default=300,
    )
    figure_fontsize = pexConfig.Field[float](
        doc="Font size for axis labels",
        default=9,
    )
    figure_height = pexConfig.Field[float](
        doc="Height of the figure. If unspecified, will be automatically set " "based on the figure width.",
        default=None,
        optional=True,
    )
    figure_width = pexConfig.Field[float](doc="Width of the figure", default=5)
    figure_histogram_height = pexConfig.Field[float](
        doc="The fraction of the bottom of the figure to be used for a"
        " histogram of metric values. If zero, it will not be shown.",
        default=0.25,
        check=lambda x: 0 <= x < 1,
    )
    key_patch = pexConfig.Field[str](doc="Key for the patch index", default="tract_patch_metric_patch_{band}")
    key_tract = pexConfig.Field[str](doc="Key for the tract index", default="tract_patch_metric_tract_{band}")
    keys_coord_ra = pexConfig.ListField[str](
        doc="Keys for the right ascension coordinates to set plot limits from "
        "if ra_min and/or ra_max are not specified.",
        default=[
            "coord_ra",
        ],
    )
    keys_coord_dec = pexConfig.ListField[str](
        doc="Keys for the declination coordinates to set plot limits from "
        "if dec_min and/or dec_max are not specified.",
        default=[
            "coord_dec",
        ],
    )
    label_tracts = pexConfig.Field[bool](doc="Whether to label tracts by ID", default=True)
    interactive = pexConfig.Field[bool](
        doc="Whether to make a figure for an interactive session by calling"
        " matplotlib.pyplot.figure. This should be False in pipelines.",
        default=False,
    )
    mag_star_bright_max = pexConfig.Field[float](
        doc="Maximum magnitude to plot bright stars",
        default=17,
    )
    parameterizedBand: bool = False
    ra_max = pexConfig.Field[float](doc="Maximum RA for the plot in degrees", optional=True)
    ra_min = pexConfig.Field[float](doc="Minimum RA for the plot in degrees", optional=True)

    def getInputSchema(self) -> KeyedDataSchema:
        keys_coord = []
        if self.ra_max is None or self.ra_min is None:
            keys_coord.append(self.keys_coord_ra)
        if self.dec_max is None or self.dec_min is None:
            keys_coord.append(self.keys_coord_dec)
        for keys in keys_coord:
            for key in keys:
                yield key, Vector
        yield self.key_patch, Vector
        yield self.key_tract, Vector

    def getOutputNames(self, config: pexConfig.Config | None = None) -> Iterable[str]:
        yield from self.metrics.keys()

    @staticmethod
    def _draw_polygon(
        ax,
        lon,
        lat,
        edgecolor="red",
        linestyle="solid",
        facecolor=None,
        **kwargs,
    ):
        """Plot a polygon from a list of lon, lat coordinates.

        This routine is a convenience wrapper around plot() and fill(), both
        of which work in geodesic (great circle) coordinates.

        Parameters
        ----------
        lon : `np.ndarray`
            Array of longitude points in polygon.
        lat : `np.ndarray`
            Array of latitude points in polygon.
        edgecolor : `str`, optional
            Color of polygon boundary.  Set to None for no boundary.
        linestyle : `str`, optional
            Line style for boundary.
        facecolor : `str`, optional
            Color of polygon face.  Set to None for no fill color.
        **kwargs : `dict`, optional
            Additional keywords passed to plot.
        """
        lines = ax.plot(
            np.concatenate((lon, [lon[0]])),
            np.concatenate((lat, [lat[0]])),
            color=edgecolor,
            linestyle=linestyle,
            **kwargs,
        )
        if facecolor is not None:
            ax.fill(lon, lat, color=facecolor, **kwargs)
        return lines

    def __call__(
        self, data: KeyedData, plotInfo: Mapping[str, Any] | None = None, **kwargs
    ) -> Mapping[str, Figure]:
        """Make the per-patch property map plot.

        Parameters
        ----------
        data : `KeyedData`
            The catalog to run the action on.
        """
        band = kwargs["band"]
        skymapInfo = kwargs["skymap"]

        key_tract = self.key_tract.format(band=band)
        key_patch = self.key_patch.format(band=band)

        cmap_patch = mpl.colormaps[self.cmap_metric]
        extremes = cmap_patch([0.04, 0.96])
        cmap_patch = cmap_patch.from_list(
            self.cmap_metric,
            cmap_patch(np.linspace(0.08, 0.92, int(round(0.84 * 256)))),
        )
        cmap_patch.set_extremes(under=extremes[0], over=extremes[1])

        tracts = np.unique(data[key_tract])

        # Note: this will not work well if the data is spanning RA=0=360.
        # Detecting this is difficult, so users should specify ra_min/max
        # and hope skyproj deals with it sensibly.
        ra_min, ra_max = self.ra_min, self.ra_max
        if ra_min is None:
            ra_min = min(np.nanmin(data[key_ra]) for key_ra in self.keys_coord_ra)
        if ra_max is None:
            ra_max = min(np.nanmax(data[key_ra]) for key_ra in self.keys_coord_ra)

        dec_min, dec_max = self.dec_min, self.dec_max
        if dec_min is None:
            dec_min = min(np.nanmin(data[key_dec]) for key_dec in self.keys_coord_dec)
        if dec_max is None:
            dec_max = min(np.nanmax(data[key_dec]) for key_dec in self.keys_coord_dec)

        deg2rad = math.pi / 180.0

        ra_med = (ra_min + ra_max) / 2.0
        # dec_med = (dec_min + dec_max)/2.
        dec_med = math.asin((np.sin(dec_min * deg2rad) + np.sin(dec_max * deg2rad)) / 2) / deg2rad

        figures = {}

        width_ratios = (1 - self.figure_colorbar_width, self.figure_colorbar_width)
        kwargs_fig = {}
        has_histogram = self.figure_histogram_height > 0
        if has_histogram:
            kwargs_fig["height_ratios"] = (1 - self.figure_histogram_height, self.figure_histogram_height)

        figure_width = self.figure_width
        figure_height = self.figure_height
        if figure_height is None:
            cos_dec_med = math.cos(dec_med * deg2rad)
            figure_height = (
                figure_width
                * (dec_max - dec_min)
                / (cos_dec_med * (ra_max - ra_min))
                / (1 - self.figure_histogram_height)
            )

        for name_fig, config_fig in self.metrics.items():
            fig = (plt.figure if self.interactive else make_figure)(
                figsize=(figure_width, figure_height),
                dpi=self.figure_dpi,
            )
            axes = fig.subplots(ncols=2, nrows=1 + has_histogram, width_ratios=width_ratios, **kwargs_fig)
            if has_histogram:
                gs = axes[1, 1].get_gridspec()
                for ax in axes[1, :]:
                    ax.remove()
                axes[1, 0] = fig.add_subplot(gs[1, :])

            fig.set_layout_engine("constrained")
            # fig.subplots_adjust(bottom=0.05, left=0.12, right=0.94, top=0.95)
            if not has_histogram:
                fig.suptitle(config_fig.description)
            fig.colorbar(
                mpl.cm.ScalarMappable(
                    norm=mpl.colors.Normalize(vmin=config_fig.vmin, vmax=config_fig.vmax, clip=False),
                    cmap=cmap_patch,
                ),
                cax=axes[0][1] if has_histogram else axes[1],
                extend="both",
                pad=0,
                fraction=self.figure_colorbar_width,
            )

            sp = skyproj.GnomonicSkyproj(
                ax=axes[0][0] if has_histogram else axes[0],
                lon_0=ra_med,
                lat_0=dec_med,
                extent=(ra_min, ra_max, dec_min, dec_max),
                pad=0,
            )
            for axis in (sp.ax.xaxis, sp.ax.yaxis):
                axis.label.set_fontsize(self.figure_fontsize)
                axis.set_tick_params(labelsize=self.figure_fontsize)
            addPlotInfo(fig, plotInfo)
            figures[name_fig] = (
                fig,
                axes,
                sp,
                config_fig.key.format(band=band),
                config_fig.vmin,
                config_fig.vmax,
                [],
            )

        for tract in tracts:
            tractInfo = skymapInfo[tract]

            in_tract = data[key_tract] == tract
            idx_tract = np.where(in_tract)[0]
            patches = data[key_patch][in_tract]

            if self.label_tracts:
                for _, _, sp, *_ in figures.values():
                    ra, dec = (coord.asDegrees() for coord in tractInfo.inner_sky_region.getCenter())
                    sp.ax.text(
                        np.clip(ra, ra_min, ra_max),
                        np.clip(dec, dec_min, dec_max),
                        tractInfo.tract_id,
                        ha="center",
                        va="center",
                        size=5,
                        c=[0.45, 0.25, 0.45],
                    )

            for idx_patch, patch in enumerate(patches):
                patchInfo = tractInfo[patch]

                vertices = patchInfo.getInnerSkyPolygon().getVertices()
                clipped = tractInfo.inner_sky_region.clipTo(
                    sphgeom.Box(sphgeom.LonLat(vertices[0]), sphgeom.LonLat(vertices[2]))
                )
                lonlats = np.array(
                    [
                        [x.asDegrees() for x in (lonlat.getA(), lonlat.getB())]
                        for lonlat in (clipped.getLon(), clipped.getLat())
                    ]
                )
                lons = np.concat((lonlats[0, :], lonlats[0, ::-1]))
                lats = np.repeat(lonlats[1, :], 2)

                for name_fig, (fig, axes, sp, quant, vmin, vmax, line_values) in figures.items():
                    value = data[quant][idx_tract[idx_patch]]

                    color = cmap_patch((value - vmin) / (vmax - vmin))
                    lines = self._draw_polygon(
                        sp.ax,
                        lons,
                        lats,
                        edgecolor="k",
                        facecolor=color,
                        linewidth=0.5,
                    )
                    line_values.append((tract, patch, lines[0], value))

        for name_fig, (fig, axes, sp, quant, vmin, vmax, line_values) in figures.items():
            patch_coordinate_entries = []
            values = []
            nonfin_count = 0
            for tract, patch, lines, value in line_values:
                transform = lines.get_transform()
                vertices = transform.transform_path(lines.get_path()).vertices
                entry = {
                    "min_x": min(vertices[:, 0]),
                    "max_x": max(vertices[:, 0]),
                    "min_y": min(vertices[:, 1]),
                    "max_y": max(vertices[:, 1]),
                    "id": f"{tract},{patch}",
                    "value": f"{value:.3}",
                }
                patch_coordinate_entries.append(entry)
                if np.isfinite(value):
                    values.append(value)
                else:
                    nonfin_count += 1
            # Other attempts to avoid having the skyproj axes and colorbar
            # labels stay within the figure boundaries don't seem to work
            # fig.tight_layout()

            if has_histogram:
                values = np.array(values)
                n_bins = max(int(np.sqrt(np.sum((values > vmin) & (values < vmax)))), 10)
                # Define bins to accumulate extreme values in outer bins
                vmargin = (vmax - vmin) / n_bins
                bins = np.concat(([vmin - vmargin], np.linspace(vmin, vmax, n_bins + 1), [vmax + vmargin]))
                axis = axes[1][0]
                try:
                    axis.hist(values, bins=bins, histtype="step")
                except ValueError:
                    message = (
                        f"Histogram with bins={','.join(f"{b}.2f" for b in bins)} from {vmin=:.2f},"
                        f" {vmargin=:.2f}, {vmax=:.2f}, {n_bins=}"
                    )
                    axis.text(
                        0.5,
                        0.5,
                        message,
                        horizontalalignment="center",
                        verticalalignment="center",
                        transform=axis.transAxes,
                    )
                    _LOG.warning(message)

                axis.set_title(self.metrics[name_fig].description, fontsize=self.figure_fontsize)

            fig.metadata = {"label": "Tract,Patch", "boxes": json.dumps(patch_coordinate_entries)}

        return {k: v[0] for k, v in figures.items()}
