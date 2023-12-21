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

import unittest

import astropy.units as u
import lsst.afw.table as afwTable
import lsst.geom
import lsst.skymap
import numpy as np
import pandas as pd
from astropy.table import Table
from lsst.analysis.tools.tasks import AstrometricCatalogMatchConfig, AstrometricCatalogMatchTask
from lsst.daf.base import PropertyList
from lsst.meas.algorithms.testUtils import MockRefcatDataId
from lsst.pipe.base import InMemoryDatasetHandle
from lsst.pipe.tasks.loadReferenceCatalog import LoadReferenceCatalogTask


class TestCatalogMatch(unittest.TestCase):
    """Test AstrometricCatalogMatchTask"""

    def setUp(self):
        config = AstrometricCatalogMatchConfig()
        config.bands = ["g", "r", "i", "z", "y"]
        self.task = AstrometricCatalogMatchTask(config=config)
        self.task.config.extraColumns.append("sourceId")
        self.task.config.referenceCatalogLoader.refObjLoader.requireProperMotion = False

        self.rng = np.random.default_rng(12345)

        self.skymap = self._make_skymap()
        self.tract = 9813

        tract = self.skymap.generateTract(self.tract)
        self.tractPoly = tract.getOuterSkyPolygon()
        self.tractBbox = self.tractPoly.getBoundingBox()

        self.nStars = 1000
        starIds = np.arange(self.nStars)
        starRas = (
            self.rng.random(self.nStars) * self.tractBbox.getWidth().asDegrees()
            + self.tractBbox.getLon().getA().asDegrees()
        )
        starDecs = (
            self.rng.random(self.nStars) * self.tractBbox.getHeight().asDegrees()
            + self.tractBbox.getLat().getA().asDegrees()
        )

        refDataId, deferredRefCat = self._make_refCat(starIds, starRas, starDecs, self.tractPoly)

        loaderTask = LoadReferenceCatalogTask(
            config=config.referenceCatalogLoader,
            dataIds=[refDataId],
            name="gaia_dr3_20230707",
            refCats=[deferredRefCat],
        )
        self.loadedRefCat = self.task._loadRefCat(loaderTask, tract)
        self.loadedRefCat["sourceId"] = np.arange(0, len(self.loadedRefCat))

        self.objectTable = self._make_objectCat(starIds, starRas, starDecs)

    def _make_skymap(self):
        """Make a testing skymap.

        Returns
        -------
        `lsst.skymap.ringsSkyMap.RingsSkyMap`
            Skymap that mimics the "hsc_rings_v1" skymap
        """
        skymap_config = lsst.skymap.ringsSkyMap.RingsSkyMapConfig()
        skymap_config.numRings = 120
        skymap_config.projection = "TAN"
        skymap_config.tractOverlap = 1.0 / 60
        skymap_config.pixelScale = 0.168
        return lsst.skymap.ringsSkyMap.RingsSkyMap(skymap_config)

    def _make_refCat(self, starIds, starRas, starDecs, poly):
        """Make a mock `deferredDatasetReference` and
        `DeferredDatasetHandle.dataId for a reference catalog.

        Parameters
        ----------
        starIds : `np.ndarray` of `int`
            Source ids for the simulated stars
        starRas : `np.ndarray` of `float`
            RAs of the simulated stars
        starDecs : `np.ndarray` of `float`
            Decs of the simulated stars
        poly : `lsst.sphgeom._sphgeom.ConvexPolygon`
            Bounding polygon containing the simulated stars

        Returns
        -------
        refDataId : `lsst.meas.algorithms.testUtils.MockRefcatDataId`
            Object that replicates the functionality of a dataId
        deferredRefCat : InMemoryDatasetHandle
            Object that replicates the functionality of a `DeferredDatasetRef`
        """
        refSchema = afwTable.SimpleTable.makeMinimalSchema()
        idKey = refSchema.addField("sourceId", type="I")
        fluxKey = refSchema.addField("phot_g_mean_flux", units="nJy", type=np.float64)
        fluxErrKey = refSchema.addField("phot_g_mean_fluxErr", units="nJy", type=np.float64)
        refCat = afwTable.SimpleCatalog(refSchema)
        ref_md = PropertyList()
        ref_md.set("REFCAT_FORMAT_VERSION", 1)
        refCat.table.setMetadata(ref_md)
        for i in range(len(starIds)):
            record = refCat.addNew()
            record.set(idKey, starIds[i])
            record.setRa(lsst.geom.Angle(starRas[i], lsst.geom.degrees))
            record.setDec(lsst.geom.Angle(starDecs[i], lsst.geom.degrees))
            record.set(fluxKey, 1)
            record.set(fluxErrKey, 0.01)
        refDataId = MockRefcatDataId(poly)
        deferredRefCat = InMemoryDatasetHandle(refCat, storageClass="SimpleCatalog", htm7="mockRefCat")
        return refDataId, deferredRefCat

    def _make_objectCat(self, starIds, starRas, starDecs):
        """Make a `pd.DataFrame` catalog with the columns needed for the
        object selector.

        Parameters
        ----------
        starIds : `np.ndarray` of `int`
            Source ids for the simulated stars
        starRas : `np.ndarray` of `float`
            RAs of the simulated stars
        starDecs : `np.ndarray` of `float`
            Decs of the simulated stars
        poly : `lsst.sphgeom._sphgeom.ConvexPolygon`
            Bounding polygon containing the simulated stars

        Returns
        -------
        sourceCat : `pd.DataFrame`
            Catalog containing the simulated stars
        """
        x = self.rng.random(self.nStars) * 4000
        y = self.rng.random(self.nStars) * 4000
        radecErr = 1.0 / (3600 * 10)  # Let random scatter be about 1/10 arcsecond
        sourceDict = {
            "sourceId": starIds,
            "coord_ra": starRas + self.rng.standard_normal(self.nStars) * radecErr,
            "coord_dec": starDecs + self.rng.standard_normal(self.nStars) * radecErr,
            "x": x,
            "y": y,
        }

        for key in [
            "r_psfFlux_flag",
            "y_extendedness_flag",
            "i_pixelFlags_saturatedCenter",
            "r_extendedness_flag",
            "y_extendedness",
            "g_extendedness_flag",
            "z_extendedness",
            "i_extendedness",
            "z_pixelFlags_saturatedCenter",
            "i_psfFlux_flag",
            "r_pixelFlags_saturatedCenter",
            "xy_flag",
            "r_extendedness",
            "y_pixelFlags_saturatedCenter",
            "i_extendedness_flag",
            "patch",
            "g_psfFlux_flag",
            "y_psfFlux_flag",
            "z_psfFlux_flag",
            "g_pixelFlags_saturatedCenter",
            "z_extendedness_flag",
            "g_extendedness",
        ]:
            sourceDict[key] = 0
        for key in ["detect_isPatchInner", "detect_isDeblendedSource"]:
            sourceDict[key] = 1
        for key in ["i_psfFlux", "g_psfFlux", "r_psfFlux", "y_psfFlux", "z_psfFlux"]:
            sourceDict[key] = 1000
        for key in ["z_psfFluxErr", "i_psfFluxErr", "r_psfFluxErr", "g_psfFluxErr", "y_psfFluxErr"]:
            sourceDict[key] = 1
        sourceCat = Table.from_pandas(pd.DataFrame(sourceDict))
        return sourceCat

    def test_setRefCat(self):
        """Test whether the objects in the reference catalog are in the
        expected footprint and that we get as many as expected
        """
        coord_ra = (self.loadedRefCat["ra"] * u.degree).to(u.radian).value
        coord_dec = (self.loadedRefCat["dec"] * u.degree).to(u.radian).value
        inFootprint = self.tractBbox.contains(coord_ra, coord_dec)
        self.assertTrue(inFootprint.all())
        self.assertEqual(len(self.loadedRefCat), self.nStars)

    def test_run(self):
        """Test whether `CatalogMatchTask` correctly associates the target and
        reference catalog.
        """
        output = self.task.run(
            catalog=self.objectTable, loadedRefCat=self.loadedRefCat, bands=self.task.config.bands
        )

        self.assertEqual(len(output.matchedCatalog), self.nStars)
        self.assertListEqual(
            list(output.matchedCatalog["sourceId_target"]),
            list(output.matchedCatalog["sourceId_ref"]),
        )


class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
