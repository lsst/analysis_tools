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
    "MakeDiffimKernelQuiverPlotVisitConfig",
    "MakeDiffimKernelQuiverPlotVisitTask",
)

import matplotlib.pyplot as plt
import numpy as np

import lsst.pex.config
import lsst.pipe.base as pipeBase
from lsst.obs.lsst import LsstCam
from lsst.pipe.base import connectionTypes


class MakeDiffimKernelQuiverPlotVisitConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("instrument", "visit"),
):
    spatiallySampledMetrics = connectionTypes.Input(
        doc="QA metrics evaluated in locations throughout the difference image.",
        name="difference_image_metrics",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
        multiple=True,
    )
    visit_summary = connectionTypes.Input(
        doc="A summary of the visit-level metadata.",
        name="preliminary_visit_summary",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
    )
    quiverKernelPlot = connectionTypes.Output(
        doc="A quiver plot of the diffim kernel centroid over the focal plane.",
        name="diffimPlots_kernelQuiver_QuiverPlot_visit",
        storageClass="Plot",
        dimensions=("instrument", "visit"),
    )


class MakeDiffimKernelQuiverPlotVisitConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=MakeDiffimKernelQuiverPlotVisitConnections
):

    pivot = lsst.pex.config.ChoiceField(
        dtype=str,
        default="middle",
        doc="Where to rotate the quiver about the sample point.",
        allowed={
            "middle": "Rotate the arrows of the quiver about their midpoint.",
            "tail": "Rotate the arrows of the quiver about their tail.",
            "tip": "Rotate the arrows of the quiver about their tip.",
        },
    )
    color = lsst.pex.config.Field(
        dtype=str,
        default="hsv",
        doc="Color map for quiver arrows",
    )
    width = lsst.pex.config.Field(
        dtype=float,
        default=0.0005,
        doc="Width of quiver arrows",
    )
    scale = lsst.pex.config.Field(
        dtype=float,
        default=1.0,
        doc="Scale factor for quiver arrows",
    )
    scaleUnits = lsst.pex.config.Field(
        dtype=str,
        default="xy",
        doc="Scale units for quiver arrows",
    )
    qKeySize = lsst.pex.config.Field(
        dtype=float,
        default=0.1,
        doc="Size of quiver key in arcsec",
    )
    detectorColors = lsst.pex.config.DictField(
        keytype=str,
        itemtype=str,
        default={"ITL": "blue", "E2V": "red"},
        doc="Colors to use for different detector types",
    )
    useDetectorColors = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Use different colors for different detector types",
    )
    directionKey = lsst.pex.config.Field(
        dtype=str,
        default="psfMatchingKernel_direction",
        doc="Column name to use for the angle of the vectors.",
    )
    lengthKey = lsst.pex.config.Field(
        dtype=str,
        default="psfMatchingKernel_length",
        doc="Column name to use for the length of the vectors.",
    )
    titleText = lsst.pex.config.Field(
        dtype=str,
        default="Diffim Kernel Quiver Plot",
        doc="Title text to display on the plot.",
    )


class MakeDiffimKernelQuiverPlotVisitTask(pipeBase.PipelineTask):

    ConfigClass = MakeDiffimKernelQuiverPlotVisitConfig
    _DefaultName = "makeDiffimKernelQuiverPlotVisitConfig"

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        """Takes a set of detector spatially sampled metrics for a visit and
        makes a quiver plot showing the spatial variation of the diffim kernel
        over the focal plane.


        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
            A butler which is specialized to operate in the context of a
            `lsst.daf.butler.Quantum`.
        inputRefs : `lsst.pipe.base.InputQuantizedConnection`
            Datastructure containing named attributes 'data and 'skymap'.
            The values of these attributes are the corresponding
            `lsst.daf.butler.DatasetRef` objects defined in the corresponding
            `PipelineTaskConnections` class.
        outputRefs : `lsst.pipe.base.OutputQuantizedConnection`
            Datastructure containing named attribute 'postageStamp'.
            The value of this attribute is the corresponding
            `lsst.daf.butler.DatasetRef` object defined in the corresponding
            `PipelineTaskConnections` class.
        """

        inputs = butlerQC.get(inputRefs)

        fig = self.run(**inputs)

        butlerQC.put(pipeBase.Struct(quiverKernelPlot=fig), outputRefs)

    def run(self, spatiallySampledMetrics, visit_summary):
        """Create a full focal plane quiver plot

        Parameters
        ----------
        spatiallySampledMetrics : `astropy.table.Table`
            Image quality metrics computed at spatially sampled locations.
        visit_summary : `lsst.afw.table.ExposureCatalog`
            Table of metadata for all exposures contained in the visit.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The finished figure.
        """
        quiverConf = {
            "pivot": self.config.pivot,
            "width": self.config.width,
            "scale": self.config.scale,
            "scale_units": self.config.scaleUnits,
        }
        if self.config.color not in plt.colormaps():
            quiverConf["color"] = self.config.color
        else:
            # Normalize directions to [0, 1] for colormap, covering -pi to pi
            norm = plt.Normalize(vmin=-np.pi, vmax=np.pi)
            cmap = plt.get_cmap(self.config.color)

        qKeyLabel = f"{self.config.qKeySize} arcsec Kernel offset direction"
        xAxisLabel = "RA (degrees)"
        yAxisLabel = "Dec (degrees)"

        if self.config.useDetectorColors:
            detectorColors = self.config.detectorColors
            camera = LsstCam().getCamera()

        fig = plt.figure(dpi=300, figsize=(20, 20))
        ax = fig.add_subplot(111)
        dec_avg = []
        for data_ref in spatiallySampledMetrics:
            detector = data_ref.dataId["detector"]

            data = data_ref.get()
            dataSelector = np.isfinite(data[self.config.directionKey]) & np.isfinite(
                data[self.config.lengthKey]
            )
            dataX = np.rad2deg(data["coord_ra"][dataSelector])
            dataY = np.rad2deg(data["coord_dec"][dataSelector])
            dec_avg.append(np.mean(data["coord_dec"][dataSelector]))
            dataA = data[self.config.directionKey][dataSelector]
            dataL = data[self.config.lengthKey][dataSelector]

            U = dataL * np.cos(dataA)
            V = dataL * np.sin(dataA)

            if "color" in quiverConf:
                q = ax.quiver(dataX, dataY, U, V, **quiverConf)
            else:
                colors = cmap(norm(dataA))
                q = ax.quiver(dataX, dataY, U, V, color=colors, **quiverConf)

            if self.config.useDetectorColors:
                detectorColor = detectorColors[camera[detector].getPhysicalType()]
            else:
                detectorColor = "blue"

            draw_detector_outlines(ax, detector, visit_summary, color=detectorColor, linewidth=0.5)
        dec_avg = np.mean(dec_avg)
        ax.quiverkey(q, 0.8, 0.95, self.config.qKeySize, qKeyLabel, labelpos="E", coordinates="axes")
        ax.invert_xaxis()
        ax.set_xlabel(xAxisLabel, fontsize=16)
        ax.set_ylabel(yAxisLabel, fontsize=16)
        ax.set_aspect(1 / np.cos(dec_avg))
        title = f"{self.config.titleText}\nvisit: {visit_summary['visit'][0]}"
        title += f" - {visit_summary.asAstropy()['band'][0]} band"
        ax.set_title(title, fontsize=16)
        plt.tight_layout()
        fig.canvas.draw()
        return fig


def draw_detector_outlines(ax, detector, visit_summary, color="red", linewidth=1):
    """Draw the outline of a detector on a plot.

    Parameters
    ----------
    ax : `matplotlib.axes.Axes`
        The plot to draw the detector outline on.
    detector : `int`
        The detector number to draw.
    visit_summary : `lsst.afw.table.ExposureCatalog`
        Table of metadata for all exposures contained in the visit.
    color : `str`, optional
        The color of the detector outline.
    linewidth : `int`, optional
        The width of the detector outline.
    """
    if detector not in visit_summary["id"]:
        return

    row = np.where(visit_summary["id"] == detector)[0][0]
    ra_corners = list(visit_summary["raCorners"][row])
    dec_corners = list(visit_summary["decCorners"][row])
    ra_center = np.mean(ra_corners)
    dec_center = np.mean(dec_corners)

    ra_corners.append(ra_corners[0])  # repeat to close loop
    dec_corners.append(dec_corners[0])  # repeat to close loop

    ax.plot(ra_corners, dec_corners, "--", color=color, alpha=0.9, linewidth=linewidth)
    ax.text(ra_center, dec_center, detector, alpha=1)
