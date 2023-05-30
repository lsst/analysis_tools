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
    "MatchedRefCoaddToolBase",
    "MatchedRefCoaddMetric",
    "MatchedRefCoaddDiffMagTool",
    "MatchedRefCoaddCModelFluxMetric",
    "MatchedRefCoaddDiffPositionTool",
    "MatchedRefCoaddPositionMetric",
)

from lsst.pex.config import ChoiceField, Field

from ..actions.plot.scatterplotWithTwoHists import ScatterPlotStatsAction, ScatterPlotWithTwoHists
from ..actions.vector.calcBinnedStats import CalcBinnedStatsAction
from ..actions.vector.mathActions import ConstantValue, DivideVector, MultiplyVector, SubtractVector
from ..actions.vector.selectors import GalaxySelector, RangeSelector, StarSelector
from ..actions.vector.vectorActions import ConvertFluxToMag, DownselectVector, LoadVector, VectorSelector
from ..interfaces import AnalysisTool, KeyedData


class MatchedRefCoaddToolBase(AnalysisTool):
    """Base tool for matched-to-reference metrics/plots on coadds.

    Metrics/plots are expected to use the reference magnitude and
    require separate star/galaxy/all source selections.

    Notes
    -----
    The tool does not use a standard coadd flag selector, because
    it is expected that the matcher has been configured to select
    appropriate candidates (and stores a match_candidate column).
    """

    def setDefaults(self):
        super().setDefaults()
        self.process.buildActions.fluxes_ref = LoadVector(vectorKey="refcat_flux_{band}")
        # TODO: Why won't vectorKey="fluxes_ref" work?
        # Does it need to be a filterAction?
        self.process.buildActions.mags_ref = ConvertFluxToMag(
            vectorKey=self.process.buildActions.fluxes_ref.vectorKey
        )

        # Select any finite extendedness (but still exclude NaNs)
        self.process.buildActions.allSelector = StarSelector(
            vectorKey="refExtendedness", extendedness_maximum=1.0
        )
        self.process.buildActions.galaxySelector = GalaxySelector(vectorKey="refExtendedness")
        self.process.buildActions.starSelector = StarSelector(vectorKey="refExtendedness")

        self.process.filterActions.xAll = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.xGalaxies = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.xStars = DownselectVector(
            vectorKey="mags_ref", selector=VectorSelector(vectorKey="starSelector")
        )


class MatchedRefCoaddMetric(MatchedRefCoaddToolBase):
    """Base tool for matched-to-reference metrics on coadds."""

    name_prefix = Field[str](default=None, doc="Prefix for metric key")
    unit = Field[str](default=None, doc="Astropy unit of y-axis values")

    _mag_low_min: int = 15
    _mag_low_max: int = 27
    _mag_interval: int = 1

    _names = ("stars", "galaxies", "all")
    _types = ("unresolved", "resolved", "all")

    def _validate(self):
        if self.name_prefix is None or self.unit is None:
            raise ValueError(
                f"{self.name_prefix=} and {self.unit=} must not be None;"
                f" did you forget to set a valid context?"
            )

    def configureMetrics(
        self,
        unit: str | None = None,
        name_prefix: str | None = None,
        name_suffix: str = "_ref_mag{minimum}",
        unit_select: str = "mag",
    ):
        """Configure metric actions and return units.

        Parameters
        ----------
        unit : `str`
            The (astropy) unit of the summary statistic metrics.
        name_prefix : `str`
            The prefix for the action (column) name.
        name_suffix : `str`
            The sufffix for the action (column) name.
        unit_select : `str`
            The (astropy) unit of the selection (x-axis) column. Default "mag".

        Returns
        -------
        units : `dict` [`str`, `str`]
            A dict of the unit (value) for each metric name (key)
        """
        unit_is_none = unit is None
        name_prefix_is_none = name_prefix is None

        if unit_is_none or name_prefix_is_none:
            if unit_is_none:
                unit = self.unit
            if name_prefix_is_none:
                name_prefix = self.name_prefix
            self._validate()
        if unit_select is None:
            unit_select = "mag"

        assert name_prefix is not None
        units = {}
        for name, name_class in zip(self._names, self._types):
            name_capital = name.capitalize()
            x_key = f"x{name_capital}"

            for minimum in range(self._mag_low_min, self._mag_low_max + 1):
                action = getattr(self.process.calculateActions, f"{name}{minimum}")
                action.selector_range = RangeSelector(
                    key=x_key,
                    minimum=minimum,
                    maximum=minimum + self._mag_interval,
                )

                action.name_prefix = name_prefix.format(name_class=name_class)
                if self.parameterizedBand:
                    action.name_prefix = f"{{band}}_{action.name_prefix}"
                action.name_suffix = name_suffix.format(minimum=minimum)

                units.update(
                    {
                        action.name_median: unit,
                        action.name_sigmaMad: unit,
                        action.name_count: "count",
                        action.name_select_minimum: unit_select,
                        action.name_select_median: unit_select,
                        action.name_select_maximum: unit_select,
                    }
                )
        return units

    def setDefaults(self):
        super().setDefaults()

        for name in self._names:
            name_capital = name.capitalize()
            for minimum in range(self._mag_low_min, self._mag_low_max + 1):
                setattr(
                    self.process.calculateActions,
                    f"{name}{minimum}",
                    CalcBinnedStatsAction(key_vector=f"y{name_capital}"),
                )

    def __call__(self, data: KeyedData, **kwargs):
        return super().__call__(data=data, **kwargs)


class MatchedRefCoaddDiffMagTool(MatchedRefCoaddToolBase):
    """Base tool for diffs between reference and measured coadd mags.

    Notes
    -----
    The default model flux is cModel.
    """

    def matchedRefDiffContext(self):
        self.process.buildActions.diff = SubtractVector(
            actionA=ConvertFluxToMag(
                vectorKey=self.process.buildActions.fluxes_meas.vectorKey, returnMillimags=True
            ),
            actionB=DivideVector(
                actionA=self.process.buildActions.mags_ref,
                # To convert to mmag
                actionB=ConstantValue(value=1e-3),
            ),
        )

    def matchedRefChiContext(self):
        self.process.buildActions.diff = DivideVector(
            actionA=SubtractVector(
                actionA=LoadVector(vectorKey=self.process.buildActions.fluxes_meas.vectorKey),
                actionB=LoadVector(vectorKey=self.process.buildActions.fluxes_ref.vectorKey),
            ),
            actionB=LoadVector(vectorKey=f"{self.process.buildActions.fluxes_meas.vectorKey}Err"),
        )

    def setDefaults(self):
        super().setDefaults()

        self.process.buildActions.fluxes_meas = LoadVector(vectorKey="{band}_cModelFlux")
        self.process.filterActions.yAll = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="starSelector")
        )


# The diamond inheritance on MatchedRefCoaddTool seems ok
class MatchedRefCoaddCModelFluxMetric(MatchedRefCoaddDiffMagTool, MatchedRefCoaddMetric):
    """Metric for diffs between reference and CModel coadd mags."""

    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.unit = "mmag"
        self.name_prefix = "photom_mag_cModelFlux_{name_class}_diff_"
        self.produce.metric.units = self.configureMetrics()

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.unit = ""
        self.name_prefix = "photom_mag_cModelFlux_{name_class}_chi_"
        self.produce.metric.units = self.configureMetrics()

    def setDefaults(self):
        super().setDefaults()


class MatchedRefCoaddDiffPositionTool(MatchedRefCoaddToolBase):
    """Base tool for diffs between reference and measured coadd astrometry."""

    scale_factor = Field[float](
        doc="The factor to multiply positions by (i.e. the pixel scale if coordinates have pixel units)",
        default=200,
    )
    variable = ChoiceField[str](
        doc="The astrometric variable to compute metrics for",
        allowed={
            "x": "x",
            "y": "y",
        },
        optional=False,
    )

    # TODO: Determine if this can be put back into setDefaults w/o this:
    # lsst.pex.config.config.FieldValidationError:
    # Field 'process.buildActions.pos_meas.vectorKey' failed validation:
    # Required value cannot be None
    def _setPos(self):
        self.process.buildActions.pos_meas = LoadVector(vectorKey=self.variable)
        self.process.buildActions.pos_ref = LoadVector(vectorKey=f"refcat_{self.variable}")

    def matchedRefDiffContext(self):
        self._setPos()
        self.process.buildActions.diff = MultiplyVector(
            actionA=ConstantValue(value=self.scale_factor),
            actionB=SubtractVector(
                actionA=self.process.buildActions.pos_meas,
                actionB=self.process.buildActions.pos_ref,
            ),
        )

    def matchedRefChiContext(self):
        self._setPos()
        self.process.buildActions.diff = DivideVector(
            actionA=SubtractVector(
                actionA=self.process.buildActions.pos_meas,
                actionB=self.process.buildActions.pos_ref,
            ),
            actionB=LoadVector(vectorKey=f"{self.process.buildActions.pos_meas.vectorKey}Err"),
        )

    def setDefaults(self):
        super().setDefaults()

        self.process.filterActions.yAll = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="allSelector")
        )
        self.process.filterActions.yGalaxies = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="galaxySelector")
        )
        self.process.filterActions.yStars = DownselectVector(
            vectorKey="diff", selector=VectorSelector(vectorKey="starSelector")
        )


class MatchedRefCoaddPositionMetric(MatchedRefCoaddDiffPositionTool, MatchedRefCoaddMetric):
    """Metric for diffs between reference and base coadd centroids."""

    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.unit = "mas"
        self.name_prefix = f"astrom_{self.variable}_{{name_class}}_diff_"
        self.produce.metric.units = self.configureMetrics()

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.unit = ""
        self.name_prefix = f"astrom_{self.variable}_{{name_class}}_chi_"
        self.produce.metric.units = self.configureMetrics()

    def setDefaults(self):
        super().setDefaults()


class MatchedRefCoaddPlot(AnalysisTool):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot = ScatterPlotWithTwoHists()

        self.produce.plot.plotTypes = ["galaxies", "stars"]
        self.produce.plot.xAxisLabel = "Reference Magnitude (mag)"


class MatchedRefCoaddCModelPlot(MatchedRefCoaddPlot):
    def setDefaults(self):
        super().setDefaults()
        self.produce.plot.magLabel = "cModel mag"

        # downselect the cModelFlux as well
        for prefix, plural in (("star", "Stars"), ("galaxy", "Galaxies")):
            for suffix in ("", "Err"):
                setattr(
                    self.process.filterActions,
                    f"{prefix}_cModelFlux{suffix}",
                    DownselectVector(
                        vectorKey=f"{{band}}_cModelFlux{suffix}",
                        selector=VectorSelector(vectorKey=f"{prefix}Selector"),
                    ),
                )

            statAction = ScatterPlotStatsAction(vectorKey=f"y{plural.capitalize()}")
            fluxType = f"{prefix}_cModelFlux"
            statAction.highSNSelector.fluxType = fluxType
            statAction.highSNSelector.threshold = 200
            statAction.lowSNSelector.fluxType = fluxType
            statAction.lowSNSelector.threshold = 10
            statAction.fluxType = fluxType
            setattr(self.process.calculateActions, plural, statAction)


class MatchedRefCoaddCModelFluxPlot(MatchedRefCoaddCModelPlot, MatchedRefCoaddDiffMagTool):
    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.produce.plot.yAxisLabel = "cModel - Reference mmag"

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.produce.plot.yAxisLabel = "chi = (cModel - Reference mag)/error"

    def setDefaults(self):
        super().setDefaults()


class MatchedRefCoaddPositionPlot(MatchedRefCoaddCModelPlot, MatchedRefCoaddDiffPositionTool):
    def matchedRefDiffContext(self):
        super().matchedRefDiffContext()
        self.produce.plot.yAxisLabel = f"{self.variable} position (pix)"

    def matchedRefChiContext(self):
        super().matchedRefChiContext()
        self.produce.plot.yAxisLabel = f"chi = (slot - Reference {self.variable} position)/error"

    def setDefaults(self):
        super().setDefaults()
