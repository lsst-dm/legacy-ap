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
        #TODO

    def tearDown(self):
        del self.sourceClusters

    def testStage(self):
        policyFile = pexPolicy.DefaultPolicyFile(
            "ap", "SourceClusterAttributesStageDictionary.paf", "policy")
        policy = pexPolicy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())

        # generate fake job identity
        jobIdentity = { "skyTileId": 0 }

        # create and populat clipboard 
        clipboard = Clipboard()
        clipboard.put(policy.get("inputKeys.jobIdentity"), jobIdentity)
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
            sourceClusterAttributes, cluster.PersistableSourceClusterVector))


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

