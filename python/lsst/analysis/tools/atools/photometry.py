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

__all__ = ("Ap12PsfSkyPlot", "PsfCModelSkyPlot", "PsfApRatio")

from ..actions.plot.histPlot import HistPanel, HistPlot
from ..actions.plot.skyPlot import SkyPlot
from ..actions.scalar.scalarActions import MeanAction, MedianAction, SigmaMadAction
from ..actions.vector import (
    CoaddPlotFlagSelector,
    DivideVector,
    ExtinctionCorrectedMagDiff,
    LoadVector,
    MagDiff,
    SnSelector,
    StarSelector,
)
from ..interfaces import AnalysisTool


class Ap12PsfSkyPlot(AnalysisTool):
    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        # Set this to an empty list to look at the band
        # the plot is being made in.
        self.prep.selectors.flagSelector.bands = []

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 300

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness"

        # TODO: Can we make these defaults somewhere?
        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"

        self.process.buildActions.zStars = ExtinctionCorrectedMagDiff()
        self.process.buildActions.zStars.magDiff.col1 = "{band}_ap12Flux"
        self.process.buildActions.zStars.magDiff.col2 = "{band}_psfFlux"

        self.process.calculateActions.median = MedianAction()
        self.process.calculateActions.median.vectorKey = "zStars"

        self.process.calculateActions.mean = MeanAction()
        self.process.calculateActions.mean.vectorKey = "zStars"

        self.process.calculateActions.sigmaMad = SigmaMadAction()
        self.process.calculateActions.sigmaMad.vectorKey = "zStars"

        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["stars"]
        self.produce.plot.plotName = "{band}_ap12-psf"
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"
        self.produce.plot.zAxisLabel = "Ap 12 - PSF [mag]"
        self.produce.plot.plotOutlines = False

        self.produce.metric.units = {"median": "mmag", "sigmaMad": "mmag", "mean": "mmag"}

        self.produce.metric.newNames = {
            "median": "{band}_ap12-psf_median",
            "mean": "{band}_ap12-psf_mean",
            "sigmaMad": "{band}_ap12-psf_sigmaMad",
        }


class PsfCModelSkyPlot(AnalysisTool):
    """Creates a skyPlot showing the difference between PSF and CModel mags"""

    def setDefaults(self):
        super().setDefaults()
        self.prep.selectors.flagSelector = CoaddPlotFlagSelector()
        self.prep.selectors.flagSelector.bands = []

        self.prep.selectors.snSelector = SnSelector()
        self.prep.selectors.snSelector.fluxType = "{band}_psfFlux"
        self.prep.selectors.snSelector.threshold = 300

        self.prep.selectors.starSelector = StarSelector()
        self.prep.selectors.starSelector.vectorKey = "{band}_extendedness"

        self.process.buildActions.xStars = LoadVector()
        self.process.buildActions.xStars.vectorKey = "coord_ra"
        self.process.buildActions.yStars = LoadVector()
        self.process.buildActions.yStars.vectorKey = "coord_dec"
        self.process.buildActions.starStatMask = SnSelector()
        self.process.buildActions.starStatMask.fluxType = "{band}_psfFlux"
        self.process.buildActions.starStatMask.threshold = 300

        self.process.buildActions.zStars = MagDiff()
        self.process.buildActions.zStars.col1 = "{band}_psfFlux"
        self.process.buildActions.zStars.col2 = "{band}_cModelFlux"

        self.process.calculateActions.median = MedianAction()
        self.process.calculateActions.median.vectorKey = "zStars"

        self.process.calculateActions.mean = MeanAction()
        self.process.calculateActions.mean.vectorKey = "zStars"

        self.process.calculateActions.sigmaMad = SigmaMadAction()
        self.process.calculateActions.sigmaMad.vectorKey = "zStars"

        self.produce.plot = SkyPlot()
        self.produce.plot.plotTypes = ["stars"]
        self.produce.plot.plotName = "{band}_psf-cModel"
        self.produce.plot.xAxisLabel = "R.A. (degrees)"
        self.produce.plot.yAxisLabel = "Dec. (degrees)"
        self.produce.plot.zAxisLabel = "PSF - cModel [mmag]"
        self.produce.plot.plotOutlines = False

        self.produce.metric.units = {"median": "mmag", "sigmaMad": "mmag", "mean": "mmag"}

        self.produce.metric.newNames = {
            "median": "{band}_psf_cModel_diff_median",
            "mean": "{band}_psf_cModel_diff_mean",
            "sigmaMad": "{band}_psf_cModel_diff_sigmaMad",
        }


class PsfApRatio(AnalysisTool):
    """Base class for plots or metrics which use PSF/Aperture Ratios."""

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.loadVectorPsf = LoadVector()
        self.process.buildActions.loadVectorAp = LoadVector()

        # assign keys for PSF and AP Flux
        self.process.buildActions.loadVectorPsf.vectorKey = "psFlux"
        self.process.buildActions.loadVectorAp.vectorKey = "apFlux"

        self.process.calculateActions.fluxRatio = DivideVector(
            actionA=self.process.buildActions.loadVectorPsf, actionB=self.process.buildActions.loadVectorAp
        )

        self.produce.plot = HistPlot()

        self.produce.panels["panel_flux"] = HistPanel()
        self.produce.panels["panel_flux"].label = "Psf/Ap Ratio"
        self.produce.panels["panel_flux"].hists = dict(fluxRatioMetric="Ratio")
