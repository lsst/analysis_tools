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
from unittest.mock import MagicMock, patch

import numpy as np
from urllib3.response import HTTPResponse

import lsst.skymap as skyMap
import lsst.utils.tests
from lsst.analysis.tools.utils import getPatchCorners, getTractCorners, http_client


class TestTractPatchUtils(lsst.utils.tests.TestCase):
    """Test to see if the tract and patch corner calculations are working.
    Includes a test case in which the tract spans the RA=360 line to ensure
    that RA wrapping is working as expected.
    """

    def setUp(self):

        # Tract 0, Patch 8 of the following SkyMap spans RA=360.
        self.tractIds = [0, 1]
        self.patchIds = [8, 0]
        skyMapConfig = skyMap.discreteSkyMap.DiscreteSkyMapConfig()
        skyMapConfig.raList = [0, 180]
        skyMapConfig.decList = [-1, 1]
        skyMapConfig.radiusList = [0.1, 0.1]
        self.skyMap = skyMap.DiscreteSkyMap(skyMapConfig)

        self.tractCorners = [
            [
                (358.8900906668926, -2.1096270745156804),
                (361.10981685029367, -2.1096270745156804),
                (361.10981685029367, 0.11009493220816949),
                (358.8900906668926, 0.11009493220816949),
            ],
            [
                (181.10981686420706, -0.11000243565052142),
                (178.89009065297807, -0.110002449565178),
                (178.88933977564471, 2.109719508502047),
                (181.1105676789934, 2.1097195571483582),
            ],
        ]
        self.patchCorners = [
            [
                (359.9999537372586, -1.7399434656613333),
                (360.37005438077574, -1.7399434656613333),
                (360.37005438077574, -1.369927758180603),
                (359.9999537372586, -1.369927758180603),
            ],
            [
                (181.10981695667596, -0.11000252814707863),
                (180.73987541079327, -0.1099561458277576),
                (180.73992023207617, 0.2600039998656897),
                (181.10988418299306, 0.25993833301062175),
            ],
        ]

    def testTractCorners(self):

        for i, tractId in enumerate(self.tractIds):
            np.testing.assert_array_almost_equal(
                getTractCorners(self.skyMap, tractId),
                self.tractCorners[i],
            )

    def testPatchCorners(self):

        for i, (tractId, patchId) in enumerate(zip(self.tractIds, self.patchIds)):
            tractInfo = self.skyMap.generateTract(tractId)
            np.testing.assert_array_almost_equal(getPatchCorners(tractInfo, patchId), self.patchCorners[i])


class TestHttpSessionAdapters(lsst.utils.tests.TestCase):
    """Tests the HTTP retry adapter.

    For common server-side HTTP failure scenarios, test the implementation of
    a retry adapter to ensure that retry scenarios are being applied,
    especially those that are non-standard or application-specific such as
    for HTTP POST, which is not retried by default.
    """

    def setUp(self):
        """For each test, create a set of HTTP status codes the client should
        encounter. These include server-side failures (5xx) and client rate
        limiting (429) before an eventual success (200).
        """
        responses = []
        self.patcher = patch("urllib3.connectionpool.HTTPConnectionPool._get_conn")
        self.mock = self.patcher.start()
        self.mock_url = "http://mock/api/resource/1"
        for code in [500, 503, 429, 200]:
            _response = MagicMock(spec=HTTPResponse())
            _response.status = code
            _response.connection = MagicMock()
            _response.headers = {}
            responses.append(_response)
        self.mock.return_value.getresponse.side_effect = responses

    def tearDown(self):
        """After each test, reset the patched mock object state."""
        self.patcher.stop()
        self.mock = None

    def testRetryServerErrorOnGet(self):
        """A GET request that fails before succeeding"""
        with http_client() as session:
            r = session.get(self.mock_url)
            assert r.ok
            assert len(r.raw.retries.history) == 3

    def testRetryServerErrorOnPost(self):
        """A POST request that fails before succeeding"""
        with http_client() as session:
            r = session.post(self.mock_url)
            assert r.ok
            assert len(r.raw.retries.history) == 3


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
