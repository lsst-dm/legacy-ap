import math
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as exceptions
import lsst.pex.policy as policy
import lsst.afw.detection as detection
import lsst.ap.cluster as cluster


def _eparams(source):
    t = source.getIxx() + source.getIyy()
    return ((source.getIxx() - source.getIyy()) / t,
            2.0 * source.getIxy() / t,
            math.sqrt(t))


class SourceClusterAttributesTestCase(unittest.TestCase):
    """Tests the source cluster attribute computation code.
    """
    def testPosition(self):
        """Tests cluster position computation.
        """
        sources = detection.SourceSet()
        for i in xrange(5):
            s = detection.Source()
            s.setRa(0.1*i)
            s.setDec(0.0)
            sources.append(s)
            if i != 2:
                s = detection.Source()
                s.setRa(0.2)
                s.setDec(0.1*(i - 2))
                sources.append(s)
        sca = cluster.SourceClusterAttributes(sources, 0, 0, 0)
        self.assertAlmostEqual(sca.getRa(), 0.2)
        self.assertAlmostEqual(sca.getDec(), 0.0)
        self.assertAlmostEqual(sca.getRaSigma(), math.sqrt(0.1/72))
        self.assertAlmostEqual(sca.getDecSigma(), math.sqrt(0.1/72))
        self.assertAlmostEqual(sca.getRaDecCov(), 0.0)

    def testTimes(self):
        """Tests observation time range computation.
        """
        sources = detection.SourceSet()
        for i in xrange(6):
            s = detection.Source()
            s.setTaiMidPoint(i)
            s.setFilterId(i / 3)
            sources.append(s)
        sca = cluster.SourceClusterAttributes(sources, 0, 0, 0)
        self.assertEqual(sca.getNumObs(), 6)
        self.assertEqual(sca.getEarliestObsTime(), 0)
        self.assertEqual(sca.getLatestObsTime(), 5)
        filters = sorted(sca.getFilterIds())
        self.assertEqual(filters, [0, 1])
        pfa = sca.getPerFilterAttributes(0)
        self.assertEqual(pfa.getFilterId(), 0)
        self.assertEqual(pfa.getNumObs(), 3)
        self.assertEqual(pfa.getEarliestObsTime(), 0)
        self.assertEqual(pfa.getLatestObsTime(), 2)
        pfa = sca.getPerFilterAttributes(1)
        self.assertEqual(pfa.getFilterId(), 1)
        self.assertEqual(pfa.getNumObs(), 3)
        self.assertEqual(pfa.getEarliestObsTime(), 3)
        self.assertEqual(pfa.getLatestObsTime(), 5)
        # Test error cehcking
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setObsTimeRange,
            pfa, 1.0, -1.0)
        self.assertRaises(
            exceptions.LsstException,
            cluster.SourceClusterAttributes.setObsTimeRange,
            sca, 1.0, -1.0)

    def testEllipticity(self):
        """Tests computation of ellipticity parameters.
        """
        sources = detection.SourceSet()
        esources = []
        for i in xrange(1, 5):
            s = detection.Source()
            s.setIxx(i)
            s.setIyy(i**2)
            s.setIxy(i)
            if i == 4:
                s.setFlagForDetection(1)
            else:
                esources.append(_eparams(s))
            sources.append(s)
        pfa = cluster.PerFilterSourceClusterAttributes(sources, 0, 1)
        self.assertEqual(pfa.getNumObs(), 4)
        self.assertEqual(pfa.getNumEllipticitySamples(), 3)
        eparams = []
        for i in xrange(3):
            eparams.append(sum(e[i] for e in esources) / len(esources))
        self.assertAlmostEqual(pfa.getE1(), eparams[0], 6)
        self.assertAlmostEqual(pfa.getE2(), eparams[1], 6)
        self.assertAlmostEqual(pfa.getRadius(), eparams[2], 6)
        uncertainties = []
        for i in xrange(3):
            u = sum((e[i] - eparams[i])**2 for e in esources)
            u = math.sqrt(u / (len(esources) * (len(esources) - 1)))
            uncertainties.append(u)
        self.assertAlmostEqual(pfa.getE1Sigma(), uncertainties[0], 6)
        self.assertAlmostEqual(pfa.getE2Sigma(), uncertainties[1], 6)
        self.assertAlmostEqual(pfa.getRadiusSigma(), uncertainties[2], 6)
        # Test with a single source
        sources.clear()
        s = detection.Source()
        s.setIxx(1.0)
        s.setIyy(2.0)
        s.setIxy(3.0)
        eparams = _eparams(s)
        pfa = cluster.PerFilterSourceClusterAttributes(s, 0, 0)
        self.assertAlmostEqual(pfa.getE1(), eparams[0], 6)
        self.assertAlmostEqual(pfa.getE2(), eparams[1], 6)
        self.assertAlmostEqual(pfa.getRadius(), eparams[2], 6)
        self.assertEqual(pfa.getE1Sigma(), None)
        self.assertEqual(pfa.getE2Sigma(), None)
        self.assertEqual(pfa.getRadiusSigma(), None)
        sources.append(s)
        pfa = cluster.PerFilterSourceClusterAttributes(sources, 0, 0)
        self.assertAlmostEqual(pfa.getE1(), eparams[0], 6)
        self.assertAlmostEqual(pfa.getE2(), eparams[1], 6)
        self.assertAlmostEqual(pfa.getRadius(), eparams[2], 6)
        self.assertEqual(pfa.getE1Sigma(), None)
        self.assertEqual(pfa.getE2Sigma(), None)
        self.assertEqual(pfa.getRadiusSigma(), None)
        # Test error checking
        for i in xrange(64):
            args = [1.0] * 6
            for j in xrange(6):
                if (i & (1 << j)) != 0:
                    args[j] = None
            if (all(a == None for a in args) or
                all(a == 1.0 for a in args) or
                all(a == 1.0 for a in args[:3]) and
                all(a == None for a in args[3:])):
                continue
            self.assertRaises(
                exceptions.LsstException,
                cluster.PerFilterSourceClusterAttributes.setEllipticity,
                pfa, *args)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumEllipticitySamples,
            pfa, -1)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumEllipticitySamples,
            pfa, cluster.PerFilterSourceClusterAttributes.NSAMPLE_MASK + 1)

    def testFlux(self):
        """Tests flux and flux uncertainty computation.
        """
        sources = detection.SourceSet()
        for i in xrange(6):
            s = detection.Source()
            s.setPsfFlux(i)
            s.setFlagForDetection(int(i == 5))
            sources.append(s)
        pfa = cluster.PerFilterSourceClusterAttributes(sources, 1, 0)
        self.assertEqual(pfa.getNumObs(), 6)
        self.assertEqual(pfa.getNumFluxSamples(), 5)
        self.assertEqual(pfa.getFlux(), 2)
        self.assertAlmostEqual(pfa.getFluxSigma(), 0.5 * math.sqrt(2.0), 6)
        # Test with a single source
        s = detection.Source()
        s.setPsfFlux(1.0)
        s.setPsfFluxErr(1.0)
        pfa = cluster.PerFilterSourceClusterAttributes(s, 0, 0)
        self.assertEqual(pfa.getFlux(), 1.0)
        self.assertEqual(pfa.getFluxSigma(), 1.0)
        self.assertEqual(pfa.getNumFluxSamples(), 1)
        sources.clear()
        sources.append(s)
        pfa = cluster.PerFilterSourceClusterAttributes(sources, 0, 0)
        self.assertEqual(pfa.getFlux(), 1.0)
        self.assertEqual(pfa.getFluxSigma(), 1.0)
        self.assertEqual(pfa.getNumFluxSamples(), 1)
        # Test with bad fluxes
        sources.clear()
        for i in xrange(6):
            s = detection.Source()
            s.setPsfFlux(i)
            s.setFlagForDetection(1)
            sources.append(s)
        pfa = cluster.PerFilterSourceClusterAttributes(sources, 1, 0)
        self.assertEqual(pfa.getNumFluxSamples(), 0)
        self.assertEqual(pfa.getFlux(), None)
        self.assertEqual(pfa.getFluxSigma(), None)
        self.assertEqual(pfa.getNumObs(), 6)
        # Test error checking
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setFlux,
                          pfa, None, 1.0)
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setFlux,
                          pfa, 1.0, -1.0)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumFluxSamples,
            pfa, -1)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumFluxSamples,
            pfa, cluster.PerFilterSourceClusterAttributes.NSAMPLE_MASK + 1)

    def testSourceClusterAttributes(self):
        """Tests features of the Python wrapper for
        lsst::ap::cluster::SourceClusterAttributes
        """
        sources = detection.SourceSet()
        eparams = []
        for i in xrange(6):
            s = detection.Source()
            s.setFilterId(i)
            s.setPsfFlux(i)
            s.setPsfFluxErr(i)
            s.setIxx(2.0*i)
            s.setIyy(i + 1.0)
            s.setIxy(0.5*i)
            eparams.append(_eparams(s))
            sources.append(s)
        sca = cluster.SourceClusterAttributes(sources, 1, 0, 0)
        filters = sorted(sca.getFilterIds())
        self.assertEqual(filters, range(6))
        filterMap = sca.getPerFilterAttributes()
        for f in filters:
            self.assertTrue(sca.hasFilter(f))
            self.assertTrue(f in filterMap)
            pfa = sca.getPerFilterAttributes(f)
            self.assertEqual(pfa.getNumObs(), 1)
            self.assertEqual(pfa.getNumFluxSamples(), 1)
            self.assertEqual(pfa.getNumEllipticitySamples(), 1)
            self.assertEqual(pfa.getFlux(), f)
            self.assertEqual(pfa.getFluxSigma(), f)
            self.assertAlmostEqual(pfa.getE1(), eparams[f][0], 6)
            self.assertAlmostEqual(pfa.getE2(), eparams[f][1], 6)
            self.assertAlmostEqual(pfa.getRadius(), eparams[f][2], 6)
            self.assertEqual(pfa.getE1Sigma(), None)
            self.assertEqual(pfa.getE2Sigma(), None)
            self.assertEqual(pfa.getRadiusSigma(), None)
        sca.clearPerFilterAttributes()
        filterMap = sca.getPerFilterAttributes()
        self.assertEqual(len(filterMap), 0)
        for f in filters:
            self.assertFalse(sca.hasFilter(f))
        sca.setPerFilterAttributes(pfa)
        self.assertTrue(sca.hasFilter(f))
        sca.removePerFilterAttributes(f)
        self.assertFalse(sca.hasFilter(f))

    def testBoostPersistence(self):
        """Tests boost persistence of
        lsst::ap::cluster::SourceClusterAttributes
        """
        pass

    def testDatabasePersistence(self):
        """Tests database persistence of
        lsst::ap::cluster::SourceClusterAttributes
        """
        pass


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = map(unittest.makeSuite,
        [SourceClusterAttributesTestCase,
         utilsTests.MemoryTestCase
        ])
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)

