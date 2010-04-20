import math
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as exceptions
import lsst.pex.policy as policy
import lsst.afw.detection as detection
import lsst.ap.cluster as cluster

class ClusterTestCase(unittest.TestCase):
    """Tests the OPTICS clustering implementation.
    """
    def testClusterEdgeCases(self):
        p = policy.Policy()
        ss = detection.SourceSet()
        s = detection.Source()
        s.setRa(0.0)
        s.setDec(0.0)
        ss.append(s)
        self.assertRaises(exceptions.LsstException, cluster.cluster, ss, p)
        p.set('epsilon', -1.0)
        p.set('minPoints', 0)
        p.set('pointsPerLeaf', 4)
        p.set('leafExtentThreshold', 7200.0)
        self.assertRaises(exceptions.LsstException, cluster.cluster, ss, p)
        p.set('epsilon', 0.0)
        p.set('minPoints', -1)
        self.assertRaises(exceptions.LsstException, cluster.cluster, ss, p)
        p.set('minPoints', 0)
        p.set('pointsPerLeaf', 0)
        self.assertRaises(exceptions.LsstException, cluster.cluster, ss, p)
        p.set('pointsPerLeaf', 4)
        c = cluster.cluster(ss, p)
        self.assertEqual(len(c), 1)
        p.set('minPoints', 2)
        c = cluster.cluster(ss, p)
        self.assertEqual(len(c), 0)
        p.set('epsilon', 3600.0) # 1-degree clustering distance
        s = detection.Source()
        s.setRa(math.radians(0.5))
        s.setDec(0.0)
        ss.append(s)
        s = detection.Source()
        s.setRa(0.0)
        s.setDec(math.radians(0.5))
        ss.append(s)
        c = cluster.cluster(ss, p)
        self.assertEqual(len(c), 1)
        p.set('minPoints', 3)
        c = cluster.cluster(ss, p)
        self.assertEqual(len(c), 0)
        p.set('minPoints', 0)
        p.set('epsilon', 1.0) # 1 arcsec clustering distance
        c = cluster.cluster(ss, p)
        self.assertEqual(len(c), 3) 
         
    def testCluster(self):
        p = policy.Policy()
        p.set('epsilon', 2000.0) # a little more than 0.5 deg
        p.set('minPoints', 2)
        p.set('pointsPerLeaf', 4)
        p.set('leafExtentThreshold', -1.0)
        ss = detection.SourceSet()
        # construct parallel streaks of sources
        for i in xrange(-2, 3):
            ra = 0.0
            for j in xrange(20):
                s = detection.Source()
                s.setObjectId(i)
                s.setRa(math.radians(ra))
                s.setDec(math.radians(i))
                ss.append(s)
                ra += 0.5
        # check that each streak results in a cluster
        clusters = cluster.cluster(ss, p)
        self.assertEqual(len(clusters), 5)
        for c in clusters:
            self.assertEqual(len(c), 20)
            i = c[0].getObjectId()
            for s in c:
                self.assertEqual(s.getObjectId(), i)


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

