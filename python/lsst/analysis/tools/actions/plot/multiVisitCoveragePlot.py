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

__all__ = ("MultiVisitCoveragePlot",)

import logging
from typing import List, Mapping, Optional, cast

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lsst.afw.cameraGeom import FOCAL_PLANE, Camera, DetectorType
from lsst.pex.config import Config, DictField, Field, ListField
from lsst.skymap import BaseSkyMap
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter

from ...interfaces import KeyedData, KeyedDataSchema, PlotAction, Scalar, Vector
from ..keyedData import KeyedDataSelectorAction
from ..vector.selectors import RangeSelector
from .plotUtils import mkColormap, plotProjectionWithBinning

log = logging.getLogger(__name__)


class MultiVisitCoveragePlot(PlotAction):
    """Plot the coverage for a set of visits."""

    plotName = Field[str](
        doc="The name for the plot.",
        optional=False,
    )
    cmapName = Field[str](
        doc="Name of the color map to be used if not using the default color-blind friendly "
        "orange/blue default (used if this is set to `None`). Any name available via "
        "`matplotlib.cm` may be used.",
        default=None,
        optional=True,
    )
    projection = Field[str](
        doc='Projection to plot. Currently only "raDec" and "focalPlane" are permitted. '
        "In either case, one point is plotted per visit/detector combination.",
        default="raDec",
    )
    nBins = Field[int](
        doc="Number of bins to use within the effective plot ranges along the spatial directions. "
        'Only used in the "raDec" projection (for the "focalPlane" projection, the binning is '
        "effectively one per detector).",
        default=25,
    )
    nPointBinThresh = Field[int](
        doc="Threshold number of data points above which binning of the data will be performed in "
        'the RA/Dec projection. If ``projection`` is "focalPlane", the per-detector nPoint '
        "threshold is nPointMinThresh/number of science detectors in the given ``camera``.",
        default=400,
    )
    unitsDict = DictField[str, str](
        doc="A dict mapping a parameter to its appropriate units (for label plotting).",
        default={
            "astromOffsetMean": "arcsec",
            "astromOffsetStd": "arcsec",
            "psfSigma": "pixels",
            "skyBg": "counts",
            "skyNoise": "counts",
            "visit": "number",
            "detector": "number",
            "zenithDistance": "deg",
            "zeroPoint": "mag",
            "ra": "deg",
            "dec": "deg",
            "xFp": "mm",
            "yFp": "mm",
            "medianE": "",
            "psfStarScaledDeltaSizeScatter": "",
        },
    )
    sortedFullBandList = ListField[str](
        doc="List of all bands that could, in principle, but do not have to, exist in data table. "
        "The sorting of the plot panels will follow this list (typically by wavelength).",
        default=["u", "g", "r", "i", "z", "y", "N921"],
    )
    bandLabelColorDict = DictField[str, str](
        doc="A dict mapping which color to use for the labels of a given band.",
        default={
            "u": "tab:purple",
            "g": "tab:blue",
            "r": "tab:green",
            "i": "gold",
            "z": "tab:orange",
            "y": "tab:red",
            "N921": "tab:pink",
        },
    )
    vectorsToLoadList = ListField[str](
        doc="List of columns to load from input table.",
        default=[
            "visitId",
            "detector",
            "band",
            "ra",
            "dec",
            "zeroPoint",
            "psfSigma",
            "skyBg",
            "astromOffsetMean",
            "psfStarDeltaE1Median",
            "psfStarDeltaE2Median",
            "psfStarScaledDeltaSizeScatter",
            "llcra",
            "lrcra",
            "ulcra",
            "urcra",
            "llcdec",
            "lrcdec",
            "ulcdec",
            "urcdec",
            "xSize",
            "ySize",
        ],
    )
    parametersToPlotList = ListField[str](
        doc="List of paramters to plot. They are plotted along rows and the columns "
        "plot these parameters for each band.",
        default=[
            "psfSigma",
            "astromOffsetMean",
            "medianE",
            "psfStarScaledDeltaSizeScatter",
            "skyBg",
            "zeroPoint",
        ],
    )
    tractsToPlotList = ListField[int](
        doc="List of tracts within which to limit the RA and Dec limits of the plot.",
        default=None,
        optional=True,
    )
    trimToTract = Field[bool](
        doc="Whether to trim the plot limits to the tract limit(s). Otherwise, plot "
        "will be trimmed to data limits (both will be expanded in the smaller range "
        "direction for an equal aspect square plot).",
        default=False,
    )
    doScatterInRaDec = Field[bool](
        doc="Whether to scatter the points in RA/Dec before plotting. This may be useful "
        "for visualization for surveys with tight dithering patterns.",
        default=False,
    )
    plotAllTractOutlines = Field[bool](
        doc="Whether to plot tract outlines for all tracts within the plot limits "
        "(regardless if they have any data in them).",
        default=False,
    )
    raDecLimitsDict = DictField[str, float](
        doc="A dict mapping the RA/Dec limits to apply to the plot. Set to `None` to use "
        "base limits on the default or the other config options. The dict must contain "
        "the keys raMin, ramax, decMin, decMax, e.g. "
        'raDecLimitsDict = {"raMin": 0, "raMax": 360, "decMin": -90, "decMax": 90}. '
        "Not compatible with ``trimToTract`` or ``tractsToPlotList`` (i.e. the latter two "
        "will be ignored if the dict is not `None`).",
        default=None,
        optional=True,
    )
    plotDetectorOutline = Field[bool](
        doc="Whether to plot a shaded outline of the detector size in the RA/Dec projection"
        "for reference. Ignored if ``projection`` is not raDec or no camera is provided "
        "in the inputs.",
        default=False,
    )

    def getInputSchema(self) -> KeyedDataSchema:
        base: list[tuple[str, type[Vector] | type[Scalar]]] = []
        for vector in self.vectorsToLoadList:
            base.append((vector, Vector))
        return base

    def __call__(self, data: KeyedData, **kwargs) -> Mapping[str, Figure] | Figure:
        self._validateInput(data, **kwargs)
        return self.makePlot(data, **kwargs)

    def _validateInput(self, data: KeyedData, **kwargs) -> None:
        """NOTE currently can only check that something is not a Scalar, not
        check that the data is consistent with Vector.
        """
        needed = self.getFormattedInputSchema(**kwargs)
        if remainder := {key.format(**kwargs) for key, _ in needed} - {
            key.format(**kwargs) for key in data.keys()
        }:
            raise ValueError(f"Task needs keys {remainder} but they were not found in input.")
        for name, typ in needed:
            isScalar = issubclass((colType := type(data[name.format(**kwargs)])), Scalar)
            if isScalar and typ != Scalar:
                raise ValueError(f"Data keyed by {name} has type {colType} but action requires type {typ}.")

        if "medianE" in self.parametersToPlotList:
            if not all(
                vector in self.vectorsToLoadList
                for vector in ["psfStarDeltaE1Median", "psfStarDeltaE2Median"]
            ):
                raise RuntimeError(
                    "If medianE is to be plotted, both psfStarDeltaE1Median and "
                    "psfStarDeltaE2Median must be included in vectorsToLoadList."
                )
        if self.raDecLimitsDict is not None:
            requiredKeys = set(["raMin", "raMax", "decMin", "decMax"])
            if requiredKeys != set(self.raDecLimitsDict.keys()):
                raise RuntimeError(
                    f"The following keys (and only these) are required in raDecLimitDict: {requiredKeys}."
                    f"The dict provided gave {set(self.raDecLimitsDict.keys())}."
                )

    def makePlot(
        self,
        data: KeyedData,
        plotInfo: Optional[Mapping[str, str]] = None,
        camera: Optional[Camera] = None,
        skymap: Optional[BaseSkyMap] = None,
        calibrateConfig: Optional[Config] = None,
        makeWarpConfig: Optional[Config] = None,
        **kwargs,
    ) -> Figure:
        """Make an Nband x Nparameter panel multi-visit coverage plot.

        The panels rows are for different bands,sorted according to the order
        in config ``sortedFullBandList``. The columns are per-parameter,
        plotted in the order given by the config ``parametersToPlotList``.
        Red "over" colors indicate thresholds in play in the data processing
        (e.g. used to select against data of sufficient quality).

        Parameters
        ----------
        data : `lsst.analysis.tools.interfaces.KeyedData`
            The key-based catalog of data to plot.
        plotInfo : `dict` [`str`], optional
            A dictionary of information about the data being plotted with (at
            least) keys:

            `"run"`
                Output run for the plots (`str`).
            `"tableName"`
                Name of the table from which results are taken (`str`).
        camera : `lsst.afw.cameraGeom.Camera`, optional
            The camera object associated with the data. This is to enable the
            conversion of to focal plane coordinates (if needed, i.e. for the
            focal plane projection version of this plot) and to obtain the
            projected (in RA/Dec) size of a typical detector (for reference in
            the raDec projection plots when requested, i.e. if the config
            ``plotDetectorOutline`` is `True`).
        skymap : `lsst.skymap.BaseSkyMap`, optional
            The sky map used for this dataset. If a specific tract[List] is
            provided, this is used to determine the appropriate RA & Dec limits
            to downselect the data to within those limits for plotting.
        calibrateConfig : `lsst.pex.config.Config`, optional
            The persisted config used in the calibration task for the given
            collection. Used to introspect threshold values used in the run.
        makeWarpConfig : `lsst.pex.config.Config`, optional
            The persisted config used in the makeWarp task for the given
            collection. Used to introspect threshold values used in the run.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The resulting figure.
        """
        mpl.rcParams["figure.dpi"] = 350
        mpl.rcParams["font.size"] = 5
        mpl.rcParams["xtick.labelsize"] = 5
        mpl.rcParams["ytick.labelsize"] = 5
        mpl.rcParams["xtick.major.size"] = 1.5
        mpl.rcParams["ytick.major.size"] = 1.5
        mpl.rcParams["xtick.major.width"] = 0.5
        mpl.rcParams["ytick.major.width"] = 0.5

        if "medianE" in self.parametersToPlotList:
            data["medianE"] = np.sqrt(
                data["psfStarDeltaE1Median"] ** 2.0 + data["psfStarDeltaE2Median"] ** 2.0
            )

        # Make sure all columns requested actually exist in ccdVisitTable.
        notFoundKeys = []
        for zKey in self.parametersToPlotList:
            if zKey not in data.keys():
                notFoundKeys.append(zKey)
        if len(notFoundKeys) > 0:
            raise RuntimeError(f"Requested column(s) {notFoundKeys} is(are) not in the data table.")

        maxMeanDistanceArcsec = calibrateConfig.astrometry.maxMeanDistanceArcsec  # type: ignore

        if makeWarpConfig is None:
            maxEllipResidual = 0.007
            maxScaledSizeScatter = 0.009
        else:
            maxEllipResidual = makeWarpConfig.select.value.maxEllipResidual  # type: ignore
            maxScaledSizeScatter = makeWarpConfig.select.value.maxScaledSizeScatter  # type: ignore

        cameraName = "" if camera is None else camera.getName()
        if self.projection == "focalPlane":
            if camera is None:
                raise RuntimeError("Must have an input camera if plotting focalPlane projection.")
            xKey = "xFp"
            yKey = "yFp"
            # Add the detector center in Focal Plane coords to the data table.
            detFpDict = {}
            for det in camera:  # type: ignore
                if det.getType() == DetectorType.SCIENCE:
                    detFpDict[det.getId()] = det.getCenter(FOCAL_PLANE)
            log.info("Number of SCIENCE detectors in {} camera = {}".format(cameraName, len(detFpDict)))
            xFpList = []
            yFpList = []
            for det in data["detector"]:  # type: ignore
                xFpList.append(detFpDict[det].getX())
                yFpList.append(detFpDict[det].getY())

            data["xFp"] = xFpList  # type: ignore
            data["yFp"] = yFpList  # type: ignore

            corners = camera[0].getCorners(FOCAL_PLANE)  # type: ignore
            xCorners, yCorners = zip(*corners)
            xScatLen = 0.4 * (np.nanmax(xCorners) - np.nanmin(xCorners))
            yScatLen = 0.4 * (np.nanmax(yCorners) - np.nanmin(yCorners))
            tractList: List[int] = []
        elif self.projection == "raDec":
            xKey = "ra"
            yKey = "dec"
            xScatLen, yScatLen = 0, 0
            # Use downselector without limits to get rid of any non-finite
            # RA/Dec entries.
            data = self._planeAreaSelector(data, xKey=xKey, yKey=yKey)
            if self.tractsToPlotList is not None and len(self.tractsToPlotList) > 0:
                tractList = self.tractsToPlotList
            else:
                ras = data["ra"]
                decs = data["dec"]
                tractList = list(set(skymap.findTractIdArray(ras, decs, degrees=True)))  # type: ignore
            log.info("List of tracts overlapping data:  {}".format(tractList))
            tractLimitsDict = self._getTractLimitsDict(skymap, tractList)
        else:
            raise ValueError("Unknown projection: {}".format(self.projection))

        perTract = True if self.tractsToPlotList is not None and self.projection == "raDec" else False

        if perTract and len(self.tractsToPlotList) > 0:
            data = self._trimDataToTracts(data, skymap)
            nData = len(cast(Vector, data[list(data.keys())[0]]))
            if nData == 0:
                raise RuntimeError(
                    "No data to plot. Did your tract selection of "
                    f"{self.tractsToPlotList} remove all data?"
                )

        if self.doScatterInRaDec:
            raRange = max(cast(Vector, data["ra"])) - min(cast(Vector, data["ra"]))
            decRange = max(cast(Vector, data["dec"])) - min(cast(Vector, data["dec"]))
            scatRad = max(0.05 * max(raRange, decRange), 0.12)  # min is of order of an LSSTCam detector
            log.info("Scattering data in RA/Dec within radius {:.3f} (deg)".format(scatRad))

        dataDf = pd.DataFrame(data)
        nDataId = len(dataDf)
        log.info("Number of detector data points: %i", nDataId)
        if nDataId == 0:
            raise RuntimeError("No data to plot. Exiting...")
        nVisit = len(set(dataDf["visitId"]))

        # Make a sorted list of the bands that exist in the data table.
        dataBandList = list(set(dataDf["band"]))
        bandList = [band for band in self.sortedFullBandList if band in dataBandList]
        missingBandList = list(set(dataBandList) - set(self.sortedFullBandList))
        if len(missingBandList) > 0:
            log.warning(
                "The band(s) {} are not included in self.sortedFullBandList. Please add them so "
                "they get sorted properly (namely by wavelength). For now, they will just be "
                "appended to the end of the list of those that could be sorted. You may also "
                "wish to give them entries in self.bandLabelColorList to specify the label color "
                "(otherwise, if not present, it defaults to teal).".format(missingBandList)
            )
            bandList.extend(missingBandList)
        log.info("Sorted list of existing bands: {}".format(bandList))

        ptSize = min(5, max(0.3, 600.0 / ((nDataId / len(bandList)) ** 0.5)))
        if self.doScatterInRaDec:
            ptSize = min(3, 0.7 * ptSize)

        nRow = len(bandList)
        nCol = len(self.parametersToPlotList)
        colMultiplier = 4 if nCol == 1 else 2.5
        fig, axes = plt.subplots(
            nRow,
            nCol,
            figsize=(int(colMultiplier * nCol), 2 * nRow),
            constrained_layout=True,
        )

        # Select reasonable limits for colormaps so they can be matched for
        # all bands.
        nPercent = max(2, int(0.02 * nDataId))
        vMinDict = {}
        vMaxDict = {}
        for zKey in self.parametersToPlotList:
            zKeySorted = dataDf[zKey].sort_values()
            zKeySorted = zKeySorted[np.isfinite(zKeySorted)]
            vMinDict[zKey] = np.nanmean(zKeySorted.head(nPercent))
            if zKey == "medianE":
                vMaxDict[zKey] = maxEllipResidual
            elif zKey == "psfStarScaledDeltaSizeScatter":
                vMaxDict[zKey] = maxScaledSizeScatter
            elif zKey == "astromOffsetMean" and self.projection != "raDec":
                vMaxDict[zKey] = min(maxMeanDistanceArcsec, 1.1 * np.nanmean(zKeySorted.tail(nPercent)))
            else:
                vMaxDict[zKey] = np.nanmean(zKeySorted.tail(nPercent))

        for iRow, band in enumerate(bandList):
            dataBand = dataDf[dataDf["band"] == band].copy()

            nDataIdBand = len(dataBand)
            nVisitBand = len(set(dataBand["visitId"]))
            if nDataIdBand < 2:
                log.warning("Fewer than 2 points to plot for {}. Skipping...".format(band))
                continue

            for zKey in self.parametersToPlotList:
                # Don't match the colorbars when it doesn't make sense to.
                if zKey in ["skyBg", "zeroPoint"]:
                    nPercent = max(2, int(0.02 * nDataIdBand))
                    zKeySorted = dataBand[zKey].sort_values()
                    zKeySorted = zKeySorted[np.isfinite(zKeySorted)]
                    vMinDict[zKey] = np.nanmean(zKeySorted.head(nPercent))
                    vMaxDict[zKey] = np.nanmean(zKeySorted.tail(nPercent))

            # Scatter the plots within the detector for focal plane plots.
            if self.doScatterInRaDec:
                scatRads = scatRad * np.sqrt(np.random.uniform(size=nDataIdBand))
                scatTheta = 2.0 * np.pi * np.random.uniform(size=nDataIdBand)
                xScatter = scatRads * np.cos(scatTheta)
                yScatter = scatRads * np.sin(scatTheta)
                xLabel = xKey + " + rand(scatter)"
                yLabel = yKey + " + rand(scatter)"
            else:
                xScatter = np.random.uniform(-xScatLen, xScatLen, len(dataBand[xKey]))
                yScatter = np.random.uniform(-yScatLen, yScatLen, len(dataBand[yKey]))
                xLabel = xKey
                yLabel = yKey
            dataBand["xScat"] = dataBand[xKey] + xScatter
            dataBand["yScat"] = dataBand[yKey] + yScatter
            # Accommodate odd number of quarter-turn rotated detectors.
            if self.projection == "focalPlane":
                for det in camera:  # type: ignore
                    if det.getOrientation().getNQuarter() % 2 != 0:
                        detId = int(det.getId())
                        xFpRot = dataBand.loc[dataBand.detector == detId, xKey]
                        yFpRot = dataBand.loc[dataBand.detector == detId, yKey]
                        xScatRot = dataBand.loc[dataBand.detector == detId, "xScat"]
                        yScatRot = dataBand.loc[dataBand.detector == detId, "yScat"]
                        dataBand.loc[dataBand.detector == detId, "xScat"] = xFpRot + (yScatRot - yFpRot)
                        dataBand.loc[dataBand.detector == detId, "yScat"] = yFpRot + (xScatRot - xFpRot)

            if band not in self.bandLabelColorDict:
                log.warning(
                    "The band {} is not included in the bandLabelColorList config. Please add it "
                    "to specify the label color (otherwise, it defaults to teal).".format(band)
                )
            color = self.bandLabelColorDict[band] if band in self.bandLabelColorDict else "teal"
            fontDict = {"fontsize": 5, "color": color}

            for iCol, zKey in enumerate(self.parametersToPlotList):
                if self.cmapName is None:
                    cmap = mkColormap(["darkOrange", "thistle", "midnightblue"])
                else:
                    cmap = mpl.cm.get_cmap(self.cmapName).copy()

                if zKey in ["medianE", "psfStarScaledDeltaSizeScatter"]:
                    cmap.set_over("red")
                elif (
                    zKey in ["astromOffsetMean"]
                    and self.projection != "raDec"
                    and vMaxDict[zKey] >= maxMeanDistanceArcsec
                ):
                    cmap.set_over("red")
                else:
                    if self.cmapName is None:
                        cmap.set_over("black")
                    else:
                        cmap.set_over("lemonchiffon")

                if zKey in ["psfSigma"]:
                    cmap.set_under("red")
                else:
                    if self.cmapName is None:
                        cmap.set_under("lemonchiffon")
                    else:
                        cmap.set_under("black")

                titleStr = "band: {} nVisit: {} nData: {}".format(band, nVisitBand, nDataIdBand)

                ax = axes[iRow, iCol] if axes.ndim > 1 else axes[max(iRow, iCol)]
                ax.set_title("{}".format(titleStr), loc="left", fontdict=fontDict, pad=2)
                ax.set_xlabel("{} ({})".format(xLabel, self.unitsDict[xKey]), labelpad=0)
                ax.set_ylabel("{} ({})".format(yLabel, self.unitsDict[yKey]), labelpad=1)
                ax.set_aspect("equal")
                ax.tick_params("x", labelrotation=45, pad=0)

                if self.projection == "focalPlane":
                    for det in camera:  # type: ignore
                        if det.getType() == DetectorType.SCIENCE:
                            corners = det.getCorners(FOCAL_PLANE)
                            xCorners, yCorners = zip(*corners)
                            xFpMin, xFpMax = min(xCorners), max(xCorners)
                            yFpMin, yFpMax = min(yCorners), max(yCorners)
                            detId = int(det.getId())
                            perDetData = dataBand[dataBand["detector"] == detId]
                            if len(perDetData) < 1:
                                log.debug("No data to plot for detector {}. Skipping...".format(detId))
                                continue
                            if sum(np.isfinite(perDetData[zKey])) < 1:
                                log.debug(
                                    "No finited data to plot for detector {}. Skipping...".format(detId)
                                )
                                continue
                            pcm = plotProjectionWithBinning(
                                ax,
                                perDetData["xScat"].values,
                                perDetData["yScat"].values,
                                perDetData[zKey].values,
                                cmap,
                                xFpMin,
                                xFpMax,
                                yFpMin,
                                yFpMax,
                                xNumBins=1,
                                fixAroundZero=False,
                                nPointBinThresh=max(1, int(self.nPointBinThresh / len(detFpDict))),
                                isSorted=False,
                                vmin=vMinDict[zKey],
                                vmax=vMaxDict[zKey],
                                scatPtSize=ptSize,
                            )
                if self.projection == "raDec":
                    raMin, raMax = min(cast(Vector, data["ra"])), max(cast(Vector, data["ra"]))
                    decMin, decMax = min(cast(Vector, data["dec"])), max(cast(Vector, data["dec"]))
                    raMin, raMax, decMin, decMax = _setLimitsToEqualRatio(raMin, raMax, decMin, decMax)
                    pcm = plotProjectionWithBinning(
                        ax,
                        dataBand["xScat"].values,
                        dataBand["yScat"].values,
                        dataBand[zKey].values,
                        cmap,
                        raMin,
                        raMax,
                        decMin,
                        decMax,
                        xNumBins=self.nBins,
                        fixAroundZero=False,
                        nPointBinThresh=self.nPointBinThresh,
                        isSorted=False,
                        vmin=vMinDict[zKey],
                        vmax=vMaxDict[zKey],
                        scatPtSize=ptSize,
                    )

                # Make sure all panels get the same axis limits and always make
                # equal axis ratio panels.
                if iRow == 0 and iCol == 0:
                    if self.raDecLimitsDict is not None and self.projection == "raDec":
                        xLimMin = self.raDecLimitsDict["raMin"]  # type: ignore
                        xLimMax = self.raDecLimitsDict["raMax"]  # type: ignore
                        yLimMin = self.raDecLimitsDict["decMin"]  # type: ignore
                        yLimMax = self.raDecLimitsDict["decMax"]  # type: ignore
                    else:
                        xLimMin, xLimMax = ax.get_xlim()
                        yLimMin, yLimMax = ax.get_ylim()
                        if self.projection == "focalPlane":
                            minDim = (
                                max(
                                    camera.getFpBBox().getWidth(),  # type: ignore
                                    camera.getFpBBox().getHeight(),  # type: ignore
                                )
                                / 2
                            )
                            xLimMin = min(-minDim, xLimMin)
                            xLimMax = max(minDim, xLimMax)
                        if self.trimToTract:
                            for tract, tractLimits in tractLimitsDict.items():
                                xLimMin = min(xLimMin, min(tractLimits["ras"]))
                                xLimMax = max(xLimMax, max(tractLimits["ras"]))
                                yLimMin = min(yLimMin, min(tractLimits["decs"]))
                                yLimMax = max(yLimMax, max(tractLimits["decs"]))
                            xDelta = xLimMax - xLimMin
                            xLimMin -= 0.04 * xDelta
                            xLimMax += 0.04 * xDelta

                        xLimMin, xLimMax, yLimMin, yLimMax = _setLimitsToEqualRatio(
                            xLimMin, xLimMax, yLimMin, yLimMax
                        )
                    limRange = xLimMax - xLimMin
                    yTickFmt = _tickFormatter(yLimMin * 10)

                    if self.plotAllTractOutlines and self.projection == "raDec":
                        allTractsList = []
                        for tractInfo in skymap:  # type: ignore
                            centerRa = tractInfo.getCtrCoord()[0].asDegrees()
                            centerDec = tractInfo.getCtrCoord()[1].asDegrees()
                            if (
                                centerRa > xLimMin
                                and centerRa < xLimMax
                                and centerDec > yLimMin
                                and centerDec < yLimMax
                            ):
                                allTractsList.append(int(tractInfo.getId()))
                        tractLimitsDict = self._getTractLimitsDict(skymap, allTractsList)

                upperHandles = []
                if self.doScatterInRaDec and self.projection == "raDec":
                    patch = mpl.patches.Circle(
                        (xLimMax - 1.5 * scatRad, yLimMax - 1.5 * scatRad),
                        radius=scatRad,
                        facecolor="gray",
                        edgecolor="None",
                        alpha=0.2,
                        label="scatter area",
                    )
                    ax.add_patch(patch)
                    upperHandles.append(patch)

                # Add a shaded area of the size of a detector for reference.
                if self.plotDetectorOutline and self.projection == "raDec":
                    if camera is None:
                        log.warning(
                            "Config plotDetectorOutline is True, but no camera was provided. "
                            "Reference detector outline will not be included in the plot."
                        )
                    else:
                        # Calculate area of polygon with known vertices.
                        x1, x2, x3, x4 = (
                            dataBand["llcra"],
                            dataBand["lrcra"],
                            dataBand["urcra"],
                            dataBand["ulcra"],
                        )
                        y1, y2, y3, y4 = (
                            dataBand["llcdec"],
                            dataBand["lrcdec"],
                            dataBand["urcdec"],
                            dataBand["ulcdec"],
                        )
                        areaDeg = (
                            np.abs(
                                (x1 * y2 - y1 * x2)
                                + (x2 * y3 - y2 * x3)
                                + (x3 * y4 - y3 * x4)
                                + (x4 * y1 - y4 * x1)
                            )
                            / 2.0
                        )
                        detScaleDeg = np.sqrt(areaDeg / (dataBand["xSize"] * dataBand["ySize"]))
                        detWidthDeg = np.nanmedian(detScaleDeg * dataBand["xSize"])
                        detHeightDeg = np.nanmedian(detScaleDeg * dataBand["ySize"])

                        patch = mpl.patches.Rectangle(
                            (xLimMax - 0.02 * limRange - detWidthDeg, yLimMin + 0.03 * limRange),
                            detWidthDeg,
                            detHeightDeg,
                            facecolor="turquoise",
                            alpha=0.3,
                            label="det size",
                        )
                        ax.add_patch(patch)
                        upperHandles.append(patch)

                if self.projection == "raDec":
                    if iRow == 0 and iCol == 0 and len(upperHandles) > 0:
                        ax.legend(
                            handles=upperHandles,
                            loc="upper left",
                            bbox_to_anchor=(0, 1.17 + 0.05 * len(upperHandles)),
                            edgecolor="black",
                            framealpha=0.4,
                            fontsize=5,
                        )

                # Overplot tract outlines if tracts were specified, but only if
                # the plot limits span fewer than the width of 5 tracts
                # (otherwise the labels will be too crowded to be useful).
                if len(tractList) > 0:
                    tractInfo = skymap[tractList[0]]  # type: ignore
                    tractWidthDeg = tractInfo.outer_sky_polygon.getBoundingBox().getWidth().asDegrees()
                    if limRange <= 5 * tractWidthDeg:
                        deltaLim = 0.05 * limRange
                        for tract, tractLimits in tractLimitsDict.items():
                            centerRa = tractLimits["center"][0]
                            centerDec = tractLimits["center"][1]
                            if (
                                centerRa > xLimMin + deltaLim
                                and centerRa < xLimMax - deltaLim
                                and centerDec > yLimMin + deltaLim
                                and centerDec < yLimMax - deltaLim
                            ):
                                ax.plot(tractLimits["ras"], tractLimits["decs"], color="dimgray", lw=0.4)
                                fontSize = 3 if limRange < 20 else 2
                                ax.annotate(
                                    str(tract),
                                    tractLimits["center"],
                                    va="center",
                                    ha="center",
                                    color="dimgray",
                                    fontsize=fontSize,
                                    annotation_clip=True,
                                    path_effects=[mpl.patheffects.withStroke(linewidth=1, foreground="w")],
                                )
                if self.projection == "raDec":
                    ax.set_xlim(xLimMax, xLimMin)
                else:
                    ax.set_xlim(xLimMin, xLimMax)
                ax.set_ylim(yLimMin, yLimMax)
                ax.yaxis.set_major_formatter(yTickFmt)

                # Get a tick formatter that will give all anticipated value
                # ranges the same length. This is so that the label padding
                # has the same effect on all colorbars.
                value = vMaxDict[zKey] if np.abs(vMaxDict[zKey]) > np.abs(vMinDict[zKey]) else vMinDict[zKey]
                if vMinDict[zKey] < 0 and np.abs(vMinDict[zKey]) >= vMaxDict[zKey] / 10:
                    value = vMinDict[zKey]
                tickFmt = _tickFormatter(value)
                cb = fig.colorbar(
                    pcm,
                    ax=ax,
                    extend="both",
                    aspect=14,
                    format=tickFmt,
                    pad=0.02,
                )

                cbLabel = zKey
                if zKey not in self.unitsDict:
                    log.warning(
                        "Data column {} does not have an entry in unitsDict config.  Units "
                        "will not be included in the colorbar text.".format(zKey)
                    )
                elif len(self.unitsDict[zKey]) > 0:
                    cbLabel = "{} ({})".format(zKey, self.unitsDict[zKey])

                cb.set_label(
                    cbLabel,
                    labelpad=-29,
                    color="black",
                    fontsize=6,
                    path_effects=[mpl.patheffects.withStroke(linewidth=1, foreground="w")],
                )

        runName = plotInfo["run"]  # type: ignore
        supTitle = "{} {} nVisit: {} nData: {}".format(runName, cameraName, nVisit, nDataId)
        if nCol == 1:
            supTitle = "{} {}\n nVisit: {} nData: {}".format(runName, cameraName, nVisit, nDataId)
        fig.suptitle(supTitle, fontsize=4 + nCol, ha="center")

        return fig

    def _getTractLimitsDict(self, skymap, tractList):
        """Return a dict containing tract limits needed for outline plotting.

        Parameters
        ----------
        skymap : `lsst.skymap.BaseSkyMap`
            The sky map used for this dataset. Used to obtain tract
            parameters.
        tractList : `list` [`int`]
            The list of tract ids (as integers) for which to determine the
            limits.

        Returns
        -------
        tractLimitsDict : `dict` [`dict`]
            A dictionary keyed on tract id. Each entry includes a `dict`
            including the tract RA corners, Dec corners, and the tract center,
            all in units of degrees. These are used for plotting the tract
            outlines.
        """
        tractLimitsDict = {}
        for tract in tractList:
            tractInfo = skymap[tract]
            tractBbox = tractInfo.outer_sky_polygon.getBoundingBox()
            tractCenter = tractBbox.getCenter()
            tractRa0 = (tractCenter[0] - tractBbox.getWidth() / 2).asDegrees()
            tractRa1 = (tractCenter[0] + tractBbox.getWidth() / 2).asDegrees()
            tractDec0 = (tractCenter[1] - tractBbox.getHeight() / 2).asDegrees()
            tractDec1 = (tractCenter[1] + tractBbox.getHeight() / 2).asDegrees()
            tractLimitsDict[tract] = {
                "ras": [tractRa0, tractRa1, tractRa1, tractRa0, tractRa0],
                "decs": [tractDec0, tractDec0, tractDec1, tractDec1, tractDec0],
                "center": [tractCenter[0].asDegrees(), tractCenter[1].asDegrees()],
            }

        return tractLimitsDict

    def _planeAreaSelector(
        self,
        data,
        xMin=np.nextafter(float("-inf"), 0),
        xMax=np.nextafter(float("inf"), 0),
        yMin=np.nextafter(float("-inf"), 0),
        yMax=np.nextafter(float("inf"), 0),
        xKey="ra",
        yKey="dec",
    ):
        """Helper function for downselecting on within an area on a plane.

        Parameters
        ----------
        data : `lsst.analysis.tools.interfaces.KeyedData`
            The key-based catalog of data to select on.
        xMin, xMax, yMin, yMax : `float`
            The min/max x/y values defining the range within which to
            down-select the data.
        xKey, yKey : `str`
            The column keys defining the "x" and "y" positions on the plane.

        Returns
        -------
        downSelectedData : `lsst.analysis.tools.interfaces.KeyedData`
            The down-selected catalog.
        """
        xSelector = RangeSelector(key=xKey, minimum=xMin, maximum=xMax)
        ySelector = RangeSelector(key=yKey, minimum=yMin, maximum=yMax)
        keyedSelector = KeyedDataSelectorAction(vectorKeys=data.keys())
        keyedSelector.selectors.xSelector = xSelector
        keyedSelector.selectors.ySelector = ySelector
        downSelectedData = keyedSelector(data)

        return downSelectedData

    def _trimDataToTracts(self, data, skymap):
        """Trim the data to limits set by tracts in self.tractsToPlotList.

        Parameters
        ----------
        data : `lsst.analysis.tools.interfaces.KeyedData`
            The key-based catalog of data to select on.
        skymap : `lsst.skymap.BaseSkyMap`
            The sky map used for this dataset. Used to obtain tract
            parameters.

        Returns
        -------
        downSelectedData : `lsst.analysis.tools.interfaces.KeyedData`
            The down-selected catalog.
        """
        tractRaMin = 1e12
        tractRaMax = -1e12
        tractDecMin = 1e12
        tractDecMax = -1e12
        for tract in self.tractsToPlotList:
            tractInfo = skymap[tract]
            tractBbox = tractInfo.outer_sky_polygon.getBoundingBox()
            tractCenter = tractBbox.getCenter()
            tractRa0 = (tractCenter[0] - tractBbox.getWidth() / 2).asDegrees()
            tractRa1 = (tractCenter[0] + tractBbox.getWidth() / 2).asDegrees()
            tractDec0 = (tractCenter[1] - tractBbox.getHeight() / 2).asDegrees()
            tractDec1 = (tractCenter[1] + tractBbox.getHeight() / 2).asDegrees()
            tractRaMin = min(tractRaMin, tractRa0)
            tractRaMax = max(tractRaMax, tractRa1)
            tractDecMin = min(tractDecMin, tractDec0)
            tractDecMax = max(tractDecMax, tractDec1)
        downSelectedData = self._planeAreaSelector(
            data,
            xMin=tractRaMin,
            xMax=tractRaMax,
            yMin=tractDecMin,
            yMax=tractDecMax,
            xKey="ra",
            yKey="dec",
        )
        return downSelectedData


def _tickFormatter(value):
    """Create a tick formatter such that all anticipated values end up with the
    same length.

    This accommodates values ranging from +/-0.0001 -> +/-99999

    Parameters
    ----------
    value : `float`
        The the value used to determine the appropriate formatting.

    Returns
    -------
    tickFmt : `matplotlib.ticker.FormatStrFormatter`
        The tick formatter to use with matplotlib functions.
    """
    if np.abs(value) >= 10000:
        tickFmt = FormatStrFormatter("%.0f")
    elif np.abs(value) >= 1000:
        tickFmt = FormatStrFormatter("%.1f")
        if value < 0:
            tickFmt = FormatStrFormatter("%.0f")
    elif np.abs(value) >= 100:
        tickFmt = FormatStrFormatter("%.2f")
        if value < 0:
            tickFmt = FormatStrFormatter("%.1f")
    elif np.abs(value) >= 10:
        tickFmt = FormatStrFormatter("%.3f")
        if value < 0:
            tickFmt = FormatStrFormatter("%.2f")
    elif np.abs(value) >= 1:
        tickFmt = FormatStrFormatter("%.4f")
        if value < 0:
            tickFmt = FormatStrFormatter("%.3f")
    else:
        tickFmt = FormatStrFormatter("%.4f")
        if value < 0:
            tickFmt = FormatStrFormatter("%.3f")
    return tickFmt


def _setLimitsToEqualRatio(xMin, xMax, yMin, yMax):
    """For a given set of x/y min/max, redefine to have equal aspect ratio.

    The limits are extended on both ends such that the central value is
    preserved.

    Parameters
    ----------
    xMin, xMax, yMin, yMax : `float`
        The min/max values of the x/y ranges for which to match in dynamic
        range while perserving the central values.

    Returns
    -------
    xMin, xMax, yMin, yMax : `float`
        The adjusted min/max values of the x/y ranges with equal aspect ratios.
    """
    xDelta = xMax - xMin
    yDelta = yMax - yMin
    deltaDiff = yDelta - xDelta
    if deltaDiff > 0:
        xMin -= 0.5 * deltaDiff
        xMax += 0.5 * deltaDiff
    elif deltaDiff < 0:
        yMin -= 0.5 * np.abs(deltaDiff)
        yMax += 0.5 * np.abs(deltaDiff)
    return xMin, xMax, yMin, yMax
