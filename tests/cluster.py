# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2012 LSST Corporation.
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
import lsst.pex.exceptions as ex
import lsst.afw.table as table
import lsst.ap.cluster as cluster
import lsst.afw.geom as afwGeom

def countClusters(clusters):
    """Count the number of source clusters containing more than one source.
    """
    return sum(map(lambda c: len(c) > 1, clusters))


class ClusterTestCase(unittest.TestCase):
    """Test the OPTICS clustering implementation.
    """
    def testClusterEdgeCases(self):
        cfg = cluster.ClusteringConfig()
        st = table.SourceTable.make(table.SourceTable.makeMinimalSchema())
        cat = table.SourceCatalog(st)
        cfg.epsilonArcsec = -1.0
        cfg.minNeighbors = 0
        cfg.pointsPerLeaf = 4
        cfg.leafExtentThresholdArcsec = 7200.0
        self.assertRaises(ex.Exception, cluster.cluster, cat, cfg.makeControl())
        cfg.epsilonArcsec = 0.0
        cfg.minNeighbors = -1
        self.assertRaises(ex.Exception, cluster.cluster, cat, cfg.makeControl())
        cfg.minNeighbors = 0
        cfg.pointsPerLeaf = 0
        self.assertRaises(ex.Exception, cluster.cluster, cat, cfg.makeControl())
        cfg.pointsPerLeaf = 4
        s = cat.addNew()
        s.setRa(0.0 * afwGeom.radians)
        s.setDec(0.0 * afwGeom.radians)
        c = cluster.cluster(cat, cfg.makeControl())
        self.assertEqual(len(c), 1)
        cfg.minNeighbors = 2
        c = cluster.cluster(cat, cfg.makeControl())
        self.assertEqual(countClusters(c), 0)
        cfg.epsilonArcsec = 3600.0 # 1-degree clustering distance
        s = cat.addNew()
        s.setRa(0.5 * afwGeom.degrees)
        s.setDec(0. * afwGeom.degrees)
        s = cat.addNew()
        s.setRa(0. * afwGeom.degrees)
        s.setDec(0.5 * afwGeom.degrees)
        c = cluster.cluster(cat, cfg.makeControl())
        self.assertEqual(countClusters(c), 1)
        cfg.minNeighbors = 3
        c = cluster.cluster(cat, cfg.makeControl())
        self.assertEqual(countClusters(c), 0)
        cfg.minNeighbors = 0
        cfg.epsilonArcsec = 1.0 # 1 arcsec clustering distance
        c = cluster.cluster(cat, cfg.makeControl())
        self.assertEqual(len(c), 3) 
         
    def testCluster(self):
        cfg = cluster.ClusteringConfig()
        st = table.SourceTable.make(table.SourceTable.makeMinimalSchema())
        cat = table.SourceCatalog(st)
        cfg.epsilonArcsec = 2000.0 # a little more than 0.5 deg
        cfg.minNeighbors = 2
        cfg.pointsPerLeaf = 8
        cfg.leafExtentThresholdArcsec = -1.0
        # construct 5 parallel streaks of sources
        for i in xrange(-2, 3):
            ra = 0.0
            for j in xrange(20):
                s = cat.addNew()
                s.setId(i)
                s.setRa(ra * afwGeom.degrees)
                s.setDec(i * afwGeom.degrees)
                ra += 0.5
        # check that each streak results in a cluster
        clusters = cluster.cluster(cat, cfg.makeControl())
        self.assertEqual(countClusters(clusters), 5)
        for c in clusters:
            # the 2 sources at the beginning and end of each streak may or
            # may not be assigned to a cluster
            self.assertTrue((len(c) >= 18 and len(c) <= 20) or len(c) == 1)
            i = c[0].getId()
            for s in c:
                self.assertEqual(s.getId(), i)


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = map(unittest.makeSuite,
        [ClusterTestCase,
         utilsTests.MemoryTestCase
        ])
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

