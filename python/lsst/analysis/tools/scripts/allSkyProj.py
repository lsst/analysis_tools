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

from typing import Any

import matplotlib
import matplotlib.cm as mplcm
import matplotlib.pyplot as plt
import numpy as np
import skyproj  # Eli's matplotlib + astronomical projections library
from lsst.daf.butler import Butler
from lsst.geom import SpherePoint, degrees
from lsst.sphgeom import Box, ConvexPolygon, LonLat


class RegionDisplay:
    def __init__(self, proj: skyproj.Skyproj):
        self.proj = proj
        self.bbox = Box()

    def add_polygon(self, region: ConvexPolygon, update_bbox: bool = True, **kwargs: Any) -> None:
        """Plot a sphgeom ConvexPolygon on a sky plot."""
        vertices = np.array(
            [
                (LonLat.longitudeOf(v).asDegrees(), LonLat.latitudeOf(v).asDegrees())
                for v in region.getVertices()
            ],
            dtype=float,
        )
        self.proj.fill(vertices[:, 0], vertices[:, 1], **kwargs)
        if update_bbox:
            self.bbox.expandTo(region.getBoundingBox())

    def add_polygon_corners(self, corners: list, **kwargs: Any) -> None:
        """Plot a sphgeom ConvexPolygon on a sky plot."""
        corners = np.array(corners)
        self.proj.fill(corners[:, 0], corners[:, 1], **kwargs)
        corner_coords = []
        for corner in corners:
            corner_coords.append(SpherePoint(corner[0], corner[1], units=degrees).getVector())
        polygon = ConvexPolygon(corner_coords)
        self.bbox.expandTo(polygon.getBoundingBox())

    def update_bounds(self) -> None:
        """Update the bounds of the plot to just fit all regions."""
        # print(self.bbox)
        lon = self.bbox.getLon()
        lat = self.bbox.getLat()
        self.proj.set_extent(
            [lon.getA().asDegrees(), lon.getB().asDegrees(), lat.getA().asDegrees(), lat.getB().asDegrees()]
        )


# collection: str = "HSC/runs/RC2/w_2023_07/DM-38042"
# butlerPath: str = "/repo/main"
# collection: str = "2.2i/runs/test-med-1/w_2023_13/DM-38781"
# butlerPath: str = "/repo/dc2"

collection: str = "u/csaunder/DM-39933/analysis_tools_baseline"
butlerPath: str = "/sdf/group/rubin/repo/main_20210215"

butler = Butler(butlerPath)
print("butler loaded")

dsTypes = butler.registry.queryDatasetTypes()
print("dsTypes loaded")
# print([dsType for dsType in dsTypes])
accumulator = {}
# dataID = {'instrument':'HSC', 'skymap':'hsc_rings_v1', 'band':'i',
# 'visit':1228}
# m = butler.get('sourceTable_visit_gaia_dr2_20200414_match_metrics',
# collections=[collection], dataId=dataID)
# print(m)
refs = butler.registry.queryDatasets(
    "sourceTable_visit_gaia_dr2_20200414_match_metrics", collections=[collection], findFirst=True
)
for ref in refs:
    if ref.dataId["band"] == "i":
        bundle = butler.get(ref)
        # print(bundle)
        # print(ref.dataId['band'])
        dataId = ref.dataId
        coords = butler.get(
            "sourceTable_visit",
            dataId=dataId,
            collections=[collection],
            parameters={"columns": ["coord_ra", "coord_dec"]},
        )
        min_ra = np.nanmin(coords["coord_ra"])
        max_ra = np.nanmax(coords["coord_ra"])
        min_dec = np.nanmin(coords["coord_dec"])
        max_dec = np.nanmax(coords["coord_dec"])
        corners = [(min_ra, min_dec), (max_ra, min_dec), (max_ra, max_dec), (min_ra, max_dec)]
        for name, measurements in bundle.items():
            for measurement in measurements:
                compoundName = name
                # print(name)
                container = accumulator.setdefault(compoundName, list())
                container.append((measurement.quantity.value, corners))
        # print(accumulator)
"""for dsType in dsTypes:
    if dsType.storageClass_name != "MetricValue":#"MetricMeasurementBundle":
        continue
    if "tract" not in dsType.dimensions:
        continue
    refs = list(butler.registry.queryDatasets(dsType, collections=[collection],
    findFirst=True).expanded())
    #print(refs)
"""
# for ref in refs:
# print(ref)
# bundle = butler.get(ref)
# tmp.append(bundle.metric_name)
# print('bundle.name: ', bundle.metric_name)
# help(bundle)
# for name, measurements in bundle.items():
# if 'matchedVisitCore' in str(bundle.metric_name):
#    #print(bundle.metric_name)
#    """for measurement in measurements:
#        compoundName = f"{dsType.name}_{name}_{measurement.metric_name.metric}
# "
#        container = accumulator.setdefault(compoundName, list())
#        container.append((measurement.quantity.value, ref.dataId.region))"""
# print(set(tmp))


for name, values in accumulator.items():
    minimum = min(v for v, _ in values)
    maximum = max(v for v, _ in values)
    norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum, clip=True)
    mapper = mplcm.ScalarMappable(norm=norm)

    fig, ax = plt.subplots(figsize=(8, 5))
    sp = skyproj.McBrydeSkyproj(ax=ax)
    regionDisplay = RegionDisplay(sp)
    for value, region in values:
        regionDisplay.add_polygon_corners(region, color=mapper.to_rgba(value))
    regionDisplay.update_bounds()
    plt.colorbar(mappable=mapper, ax=ax)
    plt.title(name)
    # plt.savefig(f"/sdf/home/n/nate2/testplots2/{name}_skyproj.png")
    plt.savefig(f"/sdf/home/m/mccann/rubin-user/DM-40203/testplots/{name}_skyproj.png")
    print(str(name) + "_skyproj.png saved")

    plt.close()
