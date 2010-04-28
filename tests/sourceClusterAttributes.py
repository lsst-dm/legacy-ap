import math
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as exceptions
import lsst.pex.policy as policy
import lsst.afw.detection as detection
import lsst.ap.cluster as cluster


def _ellipticity(source):
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
        pass
         
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
                esources.append(_ellipticity(s))
            sources.append(s)
        pfa = cluster.PerFilterSourceClusterAttributes(sources, 0, 1)
        self.assertEqual(pfa.getNumObs(), 4)
        eparams = []
        for i in xrange(3):
            eparams.append(sum(e[i] for e in esources) / len(esources))
        self.assertAlmostEqual(pfa.getE1(), eparams[0])
        self.assertAlmostEqual(pfa.getE2(), eparams[1])
        self.assertAlmostEqual(pfa.getRadius(), eparams[2])
        uncertainties = []
        for i in xrange(3):
            u = sum((e[i] - eparams[i])**2 for e in esources)
            u = math.sqrt(u / (len(esources) * (len(esources) - 1)))
            uncertainties.append(u)
        self.assertAlmostEqual(pfa.getE1Sigma(), uncertainties[0])
        self.assertAlmostEqual(pfa.getE2Sigma(), uncertainties[1])
        self.assertAlmostEqual(pfa.getRadiusSigma(), uncertainties[2])
        # Test with a single source
        sources.clear()

        # Test error checking
        pass

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
        self.assertAlmostEqual(pfa.getFluxSigma(), 0.5 * math.sqrt(2.0))
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


    def testPerFilterSourceClusterAttributes(self):
        """Tests the Python wrapper for
        lsst::ap::cluster::PerFilterSourceClusterAttributes
        """
        pass

    def testSourceClusterAttributes(self):
        """Tests the Python wrapper for
        lsst::ap::cluster::SourceClusterAttributes
        """
        pass

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

