import math
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as detection
import lsst.skypix as skypix
import lsst.ap.cluster as cluster

from lsst.pex.harness.Clipboard import Clipboard
from lsst.pex.harness.simpleStageTester import SimpleStageTester


def countClusters(clusters):
    """Counts the number of source cluster containing more than one source.
    """
    return sum(map(lambda c: len(c) > 1, clusters))


class SourceClusteringStageTestCase(unittest.TestCase):
    """Tests the lsst.ap.cluster.SourceClusteringStage pipeline stage.
    """
    def setUp(self):
        # construct 5 parallel streaks of sources 
        self.sources = detection.SourceSet()
        for i in xrange(-2, 3):
            ra = 0.0
            for j in xrange(20):
                s = detection.Source()
                s.setObjectId(i)
                s.setRa(math.radians(ra))
                s.setDec(math.radians(i))
                self.sources.append(s)
                ra += 0.5
        self.numInputSources = len(self.sources)

    def tearDown(self):
        del self.sources

    def testStage(self):
        policyFile = pexPolicy.DefaultPolicyFile(
            "ap", "SourceClusteringStageDictionary.paf", "policy")
        policy = pexPolicy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())

        # override various policy defaults
        policy.set("sourceClusteringPolicy.epsilonArcsec", 2000.0)
        policy.set("sourceClusteringPolicy.minPoints", 2)
        policy.set("sourceClusteringPolicy.leafExtentThresholdArcsec", -1.0)
        policy.set("quadSpherePolicy.resolutionPix", 3)
        policy.set("quadSpherePolicy.paddingArcsec", 0.0)

        # generate fake pipeline trigger event
        qs = skypix.QuadSpherePixelization(3, 0.0)
        event = dafBase.PropertySet()
        event.setInt("skyTileId", qs.id(1, 1, 1))

        # create and populat clipboard 
        clipboard = Clipboard()
        clipboard.put(policy.get("inputKeys.event"), event)
        clipboard.put(policy.get("inputKeys.sources"), self.sources)

        # run the stage
        stage = cluster.SourceClusteringStage(policy)
        tester = SimpleStageTester(stage)
        output = tester.runWorker(clipboard)

        # verify output
        self.assertTrue(output.contains(policy.get("outputKeys.sourceClusters")))
        self.assertTrue(output.contains(policy.get("outputKeys.sources")))
        self.assertTrue(output.contains(
            policy.get("outputKeys.sourceClusteringPolicy")))
        sourceClusters = output.get(policy.get("outputKeys.sourceClusters"))
        sources = output.get(policy.get("outputKeys.sources"))
        self.assertTrue(isinstance(sourceClusters, list))
        self.assertTrue(isinstance(sources, detection.PersistableSourceVector))
        self.assertTrue(self.numInputSources > len(sourceClusters))
        self.assertEqual(countClusters(sourceClusters), 5)
        for c in sourceClusters:
            if len(c) > 1:
                # the 2 sources at the beginning and end of each streak may or
                # may not be assigned to a cluster
                self.assertTrue(len(c) >= 18 and len(c) <= 20)
                i = c[0].getObjectId()
                for s in c:
                    self.assertEqual(s.getObjectId(), i)


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = [unittest.makeSuite(SourceClusteringStageTestCase),
              unittest.makeSuite(utilsTests.MemoryTestCase),
             ]
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

