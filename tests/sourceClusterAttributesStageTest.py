import math
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as detection
import lsst.ap.cluster as cluster

from lsst.pex.harness.Clipboard import Clipboard
from lsst.pex.harness.simpleStageTester import SimpleStageTester


class SourceClusterAttributesStageTestCase(unittest.TestCase):
    """Tests the lsst.ap.cluster.SourceClusteringStage pipeline stage.
    """
    def setUp(self):
        self.sourceClusters = []

    def tearDown(self):
        pass

    def testStage(self):
        policyFile = pexPolicy.DefaultPolicyFile(
            "ap", "SourceClusterAttributesStageDictionary.paf", "policy")
        policy = pexPolicy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())

        # generate fake pipeline trigger event
        event = dafBase.PropertySet()
        event.setInt("skyTileId", 0)

        # create and populat clipboard 
        clipboard = Clipboard()
        clipboard.put(policy.get("inputKeys.event"), event)
        clipboard.put(policy.get("inputKeys.sourceClusters"), self.sourceClusters)

        # run the stage
        stage = cluster.SourceClusterAttributesStage(policy)
        tester = SimpleStageTester(stage)
        output = tester.runWorker(clipboard)

        # verify output
        self.assertTrue(output.contains(
            policy.get("outputKeys.sourceClusterAttributes")))
        sourceClusterAttributes = output.get(
            policy.get("outputKeys.sourceClusterAttributes"))
        self.assertTrue(isinstance(
            sourceClusterAttributes, cluster.SourceClusterAttributesSet))


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = [unittest.makeSuite(SourceClusterAttributesStageTestCase),
              unittest.makeSuite(utilsTests.MemoryTestCase),
             ]
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

