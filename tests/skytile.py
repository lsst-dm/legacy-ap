# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import math
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.afw.detection as detection
import lsst.geom as geom
import lsst.skypix as skypix
import lsst.ap.utils as utils
import lsst.afw.coord as afwCoord
from lsst.afw.geom import degrees

class SkyTileTestCase(unittest.TestCase):
    """Tests the PT1 sky-tile class.
    """
    def testPrune1(self):
        coords = []
        for ra, dec in [ (0.0, 0.0),
                         (0.1, 0.1),
                         (0.0, 90.0),
                         (0.0, -90.0) ]:
            coords.append(afwCoord.IcrsCoord(ra * degrees, dec * degrees))
        res = 3
        qs = skypix.QuadSpherePixelization(res, 0.0)
        root, x, y = 1, 1, 1
        skyTileId = qs.id(root, x, y)
        skyTile = utils.PT1SkyTile(res, root, x, y, skyTileId)
        numContained = sum(skyTile.contains(s) for s in coords)
        self.assertEqual(numContained, 2)

    def testPrune2(self):
        coarseRes = 3
        subdiv = 6
        fineRes = coarseRes * subdiv
        coarseQs = skypix.QuadSpherePixelization(coarseRes, 0.0)
        fineQs = skypix.QuadSpherePixelization(fineRes, 0.0)
        for skyTileId in coarseQs:
            root, cx, cy = coarseQs.coords(skyTileId)
            skyTile = utils.PT1SkyTile(coarseRes, root, cx, cy, skyTileId)
            for cid in coarseQs:
                root2, cx2, cy2 = coarseQs.coords(cid)
                for x in xrange(cx2 * subdiv, (cx2 + 1) * subdiv):
                    for y in xrange(cy2 * subdiv, (cy2 + 1) * subdiv):
                        pixelCenter = geom.sphericalCoords(
                            fineQs.getCenter(fineQs.id(root2, x, y)))
                        s = afwCoord.IcrsCoord(pixelCenter[0] * degrees,
                                               pixelCenter[1] * degrees)
                        if root == root2 and cx == cx2 and cy == cy2:
                           self.assertEqual(skyTile.contains(s), True)
                        else:
                           self.assertEqual(skyTile.contains(s), False)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = map(unittest.makeSuite,
        [SkyTileTestCase,
         utilsTests.MemoryTestCase
        ])
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

