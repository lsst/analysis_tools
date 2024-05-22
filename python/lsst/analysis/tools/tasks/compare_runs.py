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

import argparse
import re
from collections import defaultdict

import lsst.geom
import matplotlib.pyplot as plt
import numpy as np
from lsst.daf.butler import Butler


class CompareMetrics:

    # def __init__(self, plotAstrometryMetrics=True):
    def __init__(self, butler, collections):

        self.butler = butler
        self.collections = collections
        self.bands = ["g", "r", "i", "z", "y"]
        self.linestyles = ["solid", "dashed", "dotted", "dashdot"]
        self.colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        # self.plotAstrometryMetrics = plotAstrometryMetrics

    def plotMatchedVisitMetrics(self, savePath, tracts=None):

        metricUnits = {}
        plotTracts = []
        allMetrics = defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: np.ones(len(self.bands)) * np.nan))
        )
        for collection in self.collections:
            matchedVisitMetricRefs = self.butler.registry.queryDatasets(
                "matchedVisitCore_metrics", collections=collection, findFirst=True
            )

            for ref in matchedVisitMetricRefs:
                tract = ref.dataId["tract"]
                if (tracts is not None) and (tract not in tracts):
                    continue
                if tract not in plotTracts:
                    plotTracts.append(tract)
                matchedVisitMetrics = self.butler.get(ref)
                for value in matchedVisitMetrics.values():
                    for metric in value:
                        band, metricName = re.match(r"(\w)_(.+)", metric.metric_name.metric).groups()
                        allMetrics[metricName][collection][tract][
                            self.bands.index(band)
                        ] = metric.quantity.value
                        if metricName not in metricUnits:
                            metricUnits[metricName] = metric.quantity.unit.name

        for metricName, metricValues in allMetrics.items():
            fig, ax = plt.subplots()
            for collection, collValues in metricValues.items():
                for tract, tractValues in collValues.items():
                    ax.plot(
                        tractValues,
                        marker="+",
                        linestyle=self.linestyles[self.collections.index(collection)],
                        color=self.colors[plotTracts.index(tract)],
                        label=f"{collection}, tract {tract}",
                    )

            ax.set_xticks(range(len(self.bands)))
            ax.set_xticklabels(self.bands)
            ax.set_xlabel("Band")
            ax.set_ylabel(metricUnits[metricName])
            ax.set_title(metricName)
            ax.legend(fontsize=8)

            fig.savefig(f"{savePath}/{metricName}.png")

    def plotAstromDiffMetrics(self, savePath, percentiles=np.array([16, 50, 84])):

        visitTable = self.butler.get("visitTable", instrument="HSC", collections=self.collections[0])
        skyMap = self.butler.get("skyMap", skymap="hsc_rings_v1", collections=self.collections[0])
        allMetrics = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        metricUnits = {}
        for collection in self.collections:
            astromDiffMetricRefs = set(
                self.butler.registry.queryDatasets(
                    "sourceTable_visit_gaia_dr3_20230707_match_metrics",
                    collections=collection,
                    findFirst=True,
                )
            )
            astromColorDiffMetricRefs = set(
                self.butler.registry.queryDatasets(
                    "sourceObjectTable_metrics", collections=collection, findFirst=True
                )
            )
            astromDiffMetricRefs.update(astromColorDiffMetricRefs)
            print(len(astromDiffMetricRefs))

            for ref in astromDiffMetricRefs:
                band = ref.dataId["band"]
                visitMetrics = self.butler.get(ref)
                for value in visitMetrics.values():
                    for metric in value:
                        metricName = metric.metric_name.metric
                        allMetrics[metricName][collection][band].append(metric.quantity.value)
                        if (metricName == "AB1_Dec") & (band == "r"):
                            visMatch = visitTable.loc[ref.dataId["visit"]]
                            center = lsst.geom.SpherePoint(
                                visMatch["ra"] * lsst.geom.degrees, visMatch["dec"] * lsst.geom.degrees
                            )
                            tract = skyMap.findTract(center)
                            print(ref.dataId["visit"], tract.tract_id, metric.quantity.value)
                        if metricName not in metricUnits:
                            metricUnits[metricName] = metric.quantity.unit.name

        for metricName, collections in allMetrics.items():
            fig, ax = plt.subplots()
            for collection, collectionValues in collections.items():
                metricPercentileArray = np.zeros((len(self.bands), 3))
                for band, bandValues in collectionValues.items():
                    metricPercentileArray[self.bands.index(band)] = np.percentile(bandValues, percentiles)
                if metricName == "AB1_Dec":
                    print(metricPercentileArray)
                ax.plot(
                    metricPercentileArray[:, 0],
                    linestyle="--",
                    label=f"{collection}, 16th Percentile",
                    color=self.colors[self.collections.index(collection)],
                )
                ax.plot(
                    metricPercentileArray[:, 1],
                    linestyle="-",
                    label=f"{collection}, 50th Percentile",
                    color=self.colors[self.collections.index(collection)],
                )
                ax.plot(
                    metricPercentileArray[:, 2],
                    linestyle="--",
                    label=f"{collection} 84th Percentile",
                    color=self.colors[self.collections.index(collection)],
                )

            ax.set_xticks(range(len(self.bands)))
            ax.set_xticklabels(self.bands)
            ax.set_xlabel("Band")
            ax.set_ylabel(metricUnits[metricName])
            ax.set_title(metricName)
            ax.legend(fontsize=8)

            fig.savefig(f"{savePath}/{metricName}.png")

    @classmethod
    def make_argument_parser(cls):
        """Make the argument parser for the command-line interface."""
        parser = argparse.ArgumentParser(
            description=(
                "Build a QuantumGraph that gathers and consolidates "
                "resource usage tables from existing metadata datasets."
            ),
        )
        parser.add_argument("repo", type=str, help="Path to data repository or butler configuration.")
        parser.add_argument(
            "collections",
            type=str,
            nargs="+",
            help="Collection(s)s to search for input metadata.",
        )
        parser.add_argument(
            "savePath",
            type=str,
        )
        parser.add_argument(
            "--plotMatchedVisitMetrics",
            type=bool,
            default=False,
        )
        parser.add_argument(
            "--plotAstromDiffMetrics",
            type=bool,
            default=True,
        )
        parser.add_argument(
            "--tracts",
            type=int,
            nargs="+",
        )

        return parser

    @classmethod
    def main(cls) -> None:
        """Run the command-line interface for this quantum-graph builder.

        This function provides the implementation for the
        ``build-gather-resource-usage-qg`` script.
        """
        parser = cls.make_argument_parser()
        args = parser.parse_args()

        butler = Butler(args.repo)

        comp = cls(butler, args.collections)
        if args.plotMatchedVisitMetrics:
            comp.plotMatchedVisitMetrics(args.savePath, tracts=args.tracts)

        if args.plotAstromDiffMetrics:
            comp.plotAstromDiffMetrics(args.savePath)


if __name__ == "__main__":
    CompareMetrics.main()
