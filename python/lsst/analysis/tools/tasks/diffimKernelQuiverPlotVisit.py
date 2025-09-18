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

import lsst.pex.config
import lsst.pipe.base as pipeBase
import matplotlib.pyplot as plt
import numpy as np
from lsst.obs.lsst import LsstCam
from lsst.pipe.base import connectionTypes

# from ..actions.plot.plotUtils import addPlotInfo


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
        doc="A postagestamp composite image of the tract.",
        name="diffimPlots_kernelQuiver_QuiverPlot_visit",
        storageClass="Plot",
        dimensions=("instrument", "visit"),
    )


class MakeDiffimKernelQuiverPlotVisitConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=MakeDiffimKernelQuiverPlotVisitConnections
):

    pivot = lsst.pex.config.Field(
        dtype=str,
        default="mid",
        doc="Perform diffim decorrelation to undo pixel correlation due to A&L ",
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
    scale_units = lsst.pex.config.Field(
        dtype=str,
        default="xy",
        doc="Scale units for quiver arrows",
    )
    qKeySize = lsst.pex.config.Field(
        dtype=float,
        default=0.1,
        doc="Size of quiver key in arcsec",
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

        quiverConf = {
            "pivot": self.config.pivot,
            "width": self.config.width,
            "scale": self.config.scale,
            "scale_units": self.config.scale_units,
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

        camera = LsstCam().getCamera()
        detectorColors = {"ITL": "blue", "E2V": "red"}

        fig = plt.figure(dpi=300, figsize=(20, 20))
        ax = fig.add_subplot(111)
        dec_avg = []
        for data_ref in spatiallySampledMetrics:
            detector = data_ref.dataId["detector"]
            data = data_ref.get()
            dataSelector = np.isfinite(data["psfMatchingKernel_direction"]) & np.isfinite(
                data["psfMatchingKernel_length"]
            )
            dataX = np.rad2deg(data["coord_ra"][dataSelector])
            dataY = np.rad2deg(data["coord_dec"][dataSelector])
            dec_avg.append(np.mean(data["coord_dec"][dataSelector]))
            dataA = data["psfMatchingKernel_direction"][dataSelector]
            dataL = data["psfMatchingKernel_length"][dataSelector]

            U = dataL * np.cos(dataA)
            V = dataL * np.sin(dataA)

            if "color" in quiverConf:
                q = ax.quiver(dataX, dataY, U, V, **quiverConf)
            else:
                colors = cmap(norm(dataA))
                q = ax.quiver(dataX, dataY, U, V, color=colors, **quiverConf)
            color = detectorColors[camera[detector].getPhysicalType()]
            draw_detector_outlines(detector, visit_summary, color=color, ax=ax, linewidth=0.5)
        dec_avg = np.mean(dec_avg)
        ax.quiverkey(q, 0.8, 0.95, self.config.qKeySize, qKeyLabel, labelpos="E", coordinates="axes")
        ax.invert_xaxis()
        ax.set_xlabel(xAxisLabel)
        ax.set_ylabel(yAxisLabel)
        ax.set_aspect(1 / np.cos(dec_avg))

        plt.subplots_adjust(wspace=0.0, hspace=0.0, right=0.85)

        fig.canvas.draw()
        return fig


def draw_detector_outlines(detector, visit_summary, color="red", ax=None, linewidth=1):
    row = np.where(visit_summary["id"] == detector)[0][0]
    ra_corners = list(visit_summary["raCorners"][row])
    dec_corners = list(visit_summary["decCorners"][row])
    ra_center = np.mean(ra_corners)
    dec_center = np.mean(dec_corners)

    ra_corners.append(ra_corners[0])  # repeat to close loop
    dec_corners.append(dec_corners[0])  # repeat to close loop

    ax.plot(ra_corners, dec_corners, "--", color=color, alpha=0.9, linewidth=linewidth)
    ax.text(ra_center, dec_center, detector, alpha=1)
