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
import lsst.ap.cluster as cluster

class SkyTileTestCase(unittest.TestCase):
    """Tests the PT1 sky-tile class.
    """
    def testPrune1(self):
        ss = detection.SourceSet()
        for ra, dec in [ (0.0, 0.0),
                         (0.1, 0.1),
                         (0.0, 90.0),
                         (0.0, -90.0) ]:
            s = detection.Source()
            s.setRa(math.radians(ra)) 
            s.setDec(math.radians(dec))
            ss.append(s)
        res = 3
        qs = skypix.QuadSpherePixelization(res, 0.0)
        root, x, y = 1, 1, 1
        skyTileId = qs.id(root, x, y)
        skyTile = cluster.PT1SkyTile(res, root, x, y, skyTileId)
        total = len(ss)
        skyTile.prune(ss)
        numRemaining = len(ss)
        self.assertEqual(numRemaining, 2)

    def testPrune2(self):
        coarseRes = 3
        subdiv = 6
        fineRes = coarseRes * subdiv
        coarseQs = skypix.QuadSpherePixelization(coarseRes, 0.0)
        fineQs = skypix.QuadSpherePixelization(fineRes, 0.0)
        ss = detection.SourceSet()
        for skyTileId in coarseQs:
            root, cx, cy = coarseQs.coords(skyTileId)
            skyTile = cluster.PT1SkyTile(coarseRes, root, cx, cy, skyTileId)
            expectedRemaining = 0
            for cid in coarseQs:
                root2, cx2, cy2 = coarseQs.coords(cid)
                for x in xrange(cx2 * subdiv, (cx2 + 1) * subdiv):
                    for y in xrange(cy2 * subdiv, (cy2 + 1) * subdiv):
                        s = detection.Source()
                        pixelCenter = geom.sphericalCoords(
                            fineQs.getCenter(fineQs.id(root2, x, y)))
                        s.setRa(math.radians(pixelCenter[0]))
                        s.setDec(math.radians(pixelCenter[1]))
                        if root == root2 and cx == cx2 and cy == cy2:
                           expectedRemaining += 1
                           s.setSourceId(1)
                        else:
                           s.setSourceId(0)
                        ss.append(s)
            skyTile.prune(ss)
            self.assertEqual(len(ss), expectedRemaining)
            for s in ss:
                self.assertEqual(s.getSourceId(), 1)

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

