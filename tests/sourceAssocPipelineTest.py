#! /usr/bin/env python
from itertools import izip
import math
import optparse
import os, os.path
import pdb
import unittest

import eups

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersistence
import lsst.pex.harness as pexHarness
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import lsst.geom as geom
import lsst.skypix as skypix
import lsst.ap.cluster as apCluster

from lsst.pex.harness.Clipboard import Clipboard
from lsst.pex.harness.IOStage import InputStage, OutputStage
from lsst.pex.harness.simpleStageTester import SimpleStageTester


class SourceAssocPipelineTestCase(unittest.TestCase):
    """Tests the SourceAssoc pipeline.
    """
    def setUp(self):
        self.sources = []
        self.scPolicy = None
        self.scaPolicy = None
        self.outPolicy = None
        self.jobIdentity = dafBase.PropertySet()
        # Setup filters
        afwImage.Filter.reset()
        afwImage.Filter.define(afwImage.FilterProperty("u"), 0)
        afwImage.Filter.define(afwImage.FilterProperty("g"), 1)
        afwImage.Filter.define(afwImage.FilterProperty("r"), 2)
        afwImage.Filter.define(afwImage.FilterProperty("i"), 3)
        afwImage.Filter.define(afwImage.FilterProperty("z"), 4)
        afwImage.Filter.define(afwImage.FilterProperty("y"), 5)
        # Read input sources
        inputFiles = []
        if 'TEST_SAP_INPUT_DIR' in os.environ:
            inputDir = os.environ['TEST_SAP_INPUT_DIR']
        else:
            inputDir = os.path.join(eups.productDir("ap"), "tests")
        for root, dirs, files in os.walk(inputDir):
            for f in files:
                if f.endswith('.boost'):
                    inputFiles.append(os.path.join(root, f))
        if len(inputFiles) == 0:
            self.fail("No input sources found!")
        pol = pexPolicy.Policy()
        p = dafPersistence.Persistence.getPersistence(pol)
        dp = dafBase.PropertySet()
        for f in inputFiles:
            rsl = dafPersistence.StorageList()
            loc = dafPersistence.LogicalLocation(f)
            rsl.append(p.getRetrieveStorage("BoostStorage", loc))
            persistable = p.unsafeRetrieve("PersistableSourceVector", rsl, dp)
            self.sources.append(
                afwDetection.PersistableSourceVector.swigConvert(persistable))
        # read stage policiies
        if 'TEST_SAP_SCSTAGE_POLICY' in os.environ:
            self.scPolicy = pexPolicy.Policy.createPolicy(
                os.environ['TEST_SAP_SCSTAGE_POLICY'])
        else:
            self.scPolicy = pexPolicy.Policy()
        defScFile = pexPolicy.DefaultPolicyFile(
            "ap", "SourceClusteringStageDictionary.paf", "policy")
        defScPolicy = pexPolicy.Policy.createPolicy(
            defScFile, defScFile.getRepositoryPath())
        self.scPolicy.mergeDefaults(defScPolicy)
        if 'TEST_SAP_SCASTAGE_POLICY' in os.environ:
            self.scaPolicy = pexPolicy.Policy.createPolicy(
                os.environ['TEST_SAP_SCASTAGE_POLICY'])
        if 'TEST_SAP_OUTSTAGE_POLICY' in os.environ:
            self.outPolicy = pexPolicy.Policy.createPolicy(
                os.environ['TEST_SAP_OUTSTAGE_POLICY'])
        # Compute convex hull of input sources.
        hullVerts = []
        for psv, f in izip(self.sources, inputFiles):
            verts = []
            ss = psv.getSources()
            for s in ss:
                verts.append(geom.cartesianUnitVector(
                    (math.degrees(s.getRa()), math.degrees(s.getDec()))))
            hull = geom.convexHull(verts)
            hullVerts.extend(hull.getVertices())
        hull = geom.convexHull(hullVerts)
        print 'Hull for all sources:'
        print str(hull)
        # Choose quad-sphere resolution so that the width
        # of a sky-tile is approximately the diameter of
        # the bounding circle of the source hull
        widthDeg = 2.0 * hull.getBoundingCircle().getRadius()
        qsRes = max(int(math.floor(90.0 / widthDeg)), 3)
        qs = skypix.QuadSpherePixelization(qsRes, 0.0)
        # Pick the sky-tile having maximal intersection area with the hull
        skyTileIds = qs.intersect(hull)
        maxArea = -1.0
        skyTileId = -1
        for id in skyTileIds:
            a = hull.intersect(qs.getGeometry(id)).area()
            if a > maxArea:
                maxArea = a
                skyTileId = id
        self.scPolicy.set("quadSpherePolicy.resolutionPix", qsRes)
        self.jobIdentity.set("skyTileId", skyTileId)
        print "Quad-sphere resolution: %d" % qsRes
        print "Sky-tile id: %d" % skyTileId

    def tearDown(self):
        del self.sources
        del self.scPolicy
        del self.scaPolicy
        del self.outPolicy
        del self.jobIdentity

    def testPipeline(self):
        clipboard = Clipboard()
        clipboard.put(self.scPolicy.get("inputKeys.jobIdentity"), self.jobIdentity)
        clipboard.put(self.scPolicy.get("inputKeys.sources"), self.sources)

        # Create stages
        s1 = apCluster.SourceClusteringStage(self.scPolicy)
        t1 = SimpleStageTester(s1)
        s2 = apCluster.SourceClusterAttributesStage(self.scaPolicy)
        t2 = SimpleStageTester(s2)
        if not self.outPolicy is None:
            s3 = OutputStage(self.outPolicy)
            t3 = SimpleStageTester(s3)

        # Run the pipeline with an optional output stage
        o1 = t1.runWorker(clipboard)
        o2 = t2.runWorker(o1)
        if not self.outPolicy is None:
            o3 = t3.runWorker(o2)


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = [unittest.makeSuite(SourceAssocPipelineTestCase),
              unittest.makeSuite(utilsTests.MemoryTestCase),
             ]
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    usage = """\
usage: %prog [options]

    Reads in one or more lsst.afw.detection.PersistableSourceVectors
    persisted using BoostStorage from the given files. These are used
    as test inputs for the SourceAssoc pipeline.
    """
    parser = optparse.OptionParser(usage)
    parser.add_option(
        "-i", "--input-dir", type="string", dest="inputDir",
        help="Input directory containing boost-persisted Sources")
    parser.add_option(
        "-o", "--out-stage-policy", type="string", dest="outStagePolicy",
        help="Output stage policy file allowing persistence of test results")
    parser.add_option(
        "-c", "--sc-stage-policy", type="string", dest="scStagePolicy",
        help="Policy file for the SourceClusteringStage")
    parser.add_option(
        "-a", "--sca-stage-policy", type="string", dest="scaStagePolicy",
        help="Policy file for the SourceClusterAttributesStage")
    (options, args) = parser.parse_args()
    if options.inputDir:
        os.environ['TEST_SAP_INPUT_DIR'] = options.inputDir
    if options.outStagePolicy:
        os.environ['TEST_SAP_OUTSTAGE_POLICY'] = options.outStagePolicy
    if options.scStagePolicy:
        os.environ['TEST_SAP_SCSTAGE_POLICY'] = options.scStagePolicy
    if options.scaStagePolicy:
        os.environ['TEST_SAP_SCASTAGE_POLICY'] = options.scaStagePolicy
    run(True)

