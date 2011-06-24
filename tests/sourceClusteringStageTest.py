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
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.skypix as skypix
import lsst.ap.cluster as cluster
import lsst.ap.match as match

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
        self.sources = afwDetection.SourceSet()
        self.exposures = []
        ps = dafBase.PropertySet()
        ps.setLong("scienceCcdExposureId", 0L)
        ps.setString("FILTER", "u")
        ps.setString("TIME-MID", "2000-01-01T11:59:28.000000000Z")
        ps.setDouble("EXPTIME", 10.0)
        ps.setString("RADESYS", "FK5")
        ps.setDouble("EQUINOX", 2000.0)
        ps.setString("CTYPE1", "RA---TAN")
        ps.setString("CTYPE2", "DEC--TAN")
        ps.setString("CUNIT1", "deg")
        ps.setString("CUNIT2", "deg")
        ps.setInt("NAXIS1", 1000)
        ps.setInt("NAXIS2", 1000)
        ps.setDouble("CRPIX1", 500.5)
        ps.setDouble("CRPIX2", 500.5)
        ps.setDouble("CRVAL1", 0.0)
        ps.setDouble("CRVAL2", 0.0)
        ps.setDouble("CD1_1", 1.0/3600.0)
        ps.setDouble("CD1_2", 0.0)
        ps.setDouble("CD2_1", 0.0)
        ps.setDouble("CD2_2", 1.0/3600.0)
        self.exposures.append(ps)
        for i in xrange(-2, 3):
            ra = 0.0
            for j in xrange(20):
                s = afwDetection.Source()
                s.setObjectId(i)
                s.setAmpExposureId(0)
                s.setRa(math.radians(ra))
                s.setDec(math.radians(i))
                s.setXAstrom(ra)
                s.setYAstrom(i)
                self.sources.append(s)
                ra += 0.5
        self.numInputSources = len(self.sources)

    def tearDown(self):
        del self.sources
        del self.exposures

    def testStage(self):
        policyFile = pexPolicy.DefaultPolicyFile(
            "ap", "SourceClusteringStageDictionary.paf", "policy")
        policy = pexPolicy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())

        # override various policy defaults
        policy.set("sourceClusteringPolicy.epsilonArcsec", 2000.0)
        policy.set("sourceClusteringPolicy.minNeighbors", 2)
        policy.set("sourceClusteringPolicy.leafExtentThresholdArcsec", -1.0)
        policy.set("quadSpherePolicy.resolutionPix", 3)
        policy.set("quadSpherePolicy.paddingArcsec", 0.0)
        policy.set("debug.createGoodSourceHistogram", True)
        policy.set("debug.sourceHistogramResolution", 100)

        # generate fake job identity
        qs = skypix.QuadSpherePixelization(3, 0.0)
        jobIdentity = dafBase.PropertySet()
        jobIdentity.setInt("skyTileId", qs.id(1, 1, 1))

        # create and populat clipboard 
        clipboard = Clipboard()
        clipboard.put(policy.get("inputKeys.jobIdentity"), jobIdentity)
        clipboard.put(policy.get("inputKeys.sources"),
                      afwDetection.PersistableSourceVector(self.sources))
        clipboard.put(policy.get("inputKeys.exposures"),
                      self.exposures)

        # run the stage
        stage = cluster.SourceClusteringStage(policy)
        tester = SimpleStageTester(stage)
        output = tester.runWorker(clipboard)

        # verify output
        self.assertTrue(output.contains(
            policy.get("outputKeys.sourceClusters")))
        self.assertTrue(output.contains(policy.get("outputKeys.sources")))
        self.assertTrue(output.contains(
            policy.get("outputKeys.goodSourceHistogram")))
        self.assertTrue(output.contains(
            policy.get("outputKeys.sourceClusteringPolicy")))
        sourceClusters = output.get(policy.get("outputKeys.sourceClusters"))
        sources = output.get(policy.get("outputKeys.sources"))
        sourceClusteringPolicy = output.get(
            policy.get("outputKeys.sourceClusteringPolicy"))
        goodSourceHistogram = output.get(
            policy.get("outputKeys.goodSourceHistogram"))
        self.assertTrue(isinstance(sourceClusters, list))
        self.assertTrue(
            isinstance(sources, afwDetection.PersistableSourceVector))
        self.assertTrue(isinstance(sourceClusteringPolicy, pexPolicy.Policy))
        self.assertTrue(
            isinstance(goodSourceHistogram, afwImage.DecoratedImageU))
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

