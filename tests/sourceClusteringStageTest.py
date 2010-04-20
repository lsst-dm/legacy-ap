import math
import pdb

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.pex.harness as harness
import lsst.afw.detection as detection
import lsst.skypix as skypix
import lsst.ap.cluster as cluster

from harness.simpleStageTester import SimpleStageTester


class SourceClusteringStageTestCase(unittest.TestCase):
    """Tests the lsst.ap.cluster.SourceClusteringStage pipeline stage.
    """
    def setUp(self):
        self.sources = detection.SourceSet()
        self.numInputSources = len(self.sources)
        for i in xrange(-2, 3):
            ra = 0.0
            for j in xrange(20):
                s = detection.Source()
                s.setObjectId(i)
                s.setRa(math.radians(ra))
                s.setDec(math.radians(i))
                self.sources.append(s)
                ra += 0.5

    def tearDown(self):
        for key in self.__dict__.keys():
            del self.__dict__[key]

    def testSourceClusteringStage(self):
        policyFile = pexPolicy.DefaultPolicyFile(
            "ap", "SourceClusteringStageDictionary.paf", "policy")
        policy = pexPolicy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())
        # override various policy defaults
        policy.set("epsilonArcsec", 2000.0)
        policy.set("minPoints", 2)
        policy.set("leafExtentThresholdArcsec", -1.0)
        policy.set("resolutionPix", 3)
        policy.set("paddingArcsec", 0.0)

        # generate fake pipeline trigger event
        qs = skypix.QuadSpherePixelization(3, 0.0)
        event = dafBase.PropertySet()
        event.setInt("skyTileId", qs.id(1, 1, 1))

        # create and populat clipboard 
        clipboard = harness.Clipboard()
        clipboard.put(policy.get("inputKeys.event"), event)
        clipboard.put(policy.get("inputKeys.sources"), self.sources)

        # run the stage
        stage = cluster.SourceClusteringStage(policy)
        tester = SimpleStageTester(stage)
        output = tester.runWorker(clipboard)

        # verify output
        self.assertTrue(output.contains(policy.get("outputKeys.skyTile")))
        self.assertTrue(output.contains(policy.get("outputKeys.sourceClusters")))
        self.assertTrue(output.contains(policy.get("outputKeys.sources")))
        skyTile = output.get(policy.get("outputKeys.skyTile"))
        sourceClusters = output.get(policy.get("outputKeys.sourceClusters"))
        sources = output.get(policy.get("outputKeys.sources"))
        self.assertTrue(isinstance(skyTile, cluster.PT1SkyTile))
        self.assertTrue(isinstance(sourceClusters, list))
        self.assertTrue(isinstance(sources, detection.PersistableSourceVector))
        self.assertTrue(self.numInputSources > len(sourceClusters))
        self.assertEqual(len(sourceClusters), 5)
        for c in sourceClusters:
            self.assertEqual(len(c), 20)
            i = c[0].getObjectId()
            for s in c:
                self.assertEqual(s.getObjectId(), i)

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = map(unittest.makeSuite,
        [SourceClusteringStageTestCase,
         utilsTests.MemoryTestCase
        ])
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

