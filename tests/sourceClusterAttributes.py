from itertools import izip
import math
import os
import pdb
import platform
import tempfile
import unittest

import lsst.utils.tests as utilsTests
import lsst.pex.exceptions as exceptions
import lsst.daf.base as base
import lsst.daf.persistence as persistence
import lsst.pex.policy as policy
import lsst.afw.detection as detection
import lsst.afw.image as image
import lsst.ap.cluster as cluster


def _eparams(source):
    t = source.getIxx() + source.getIyy()
    return ((source.getIxx() - source.getIyy()) / t,
            2.0 * source.getIxy() / t,
            math.sqrt(t))


class SourceClusterAttributesTestCase(unittest.TestCase):
    """Tests the source cluster attribute computation code.
    """
    def setUp(self):
        self.clusters = cluster.SourceClusterVector()
        for i in xrange(20):
            sca = cluster.SourceClusterAttributes()
            sca.setClusterId(i)
            sca.setFlags(i)
            sca.setNumObs(0)
            sca.setObsTimeRange(float(i), float(i + 1))
            if i % 3 == 0:
                sca.setPosition(0.2*i, 0.1*(i - 10), 0.01*i, 0.05*i, 0.03*i)
            elif i % 3 == 1:
                sca.setPosition(0.2*i, 0.1*(i - 10), 0.01*i, 0.05*i, None)
            else:
                sca.setPosition(0.2*i, 0.1*(i - 10), None, None, None)
            # add in per-filter attributes
            for j in xrange(i % 7):
                pfa = cluster.PerFilterSourceClusterAttributes()
                pfa.setFilterId(j)
                pfa.setNumObs(j + 1)
                sca.setNumObs(sca.getNumObs() + pfa.getNumObs())
                pfa.setNumFluxSamples(max(0, j - 1))
                pfa.setNumFluxSamples(max(0, j - 2))
                if i % 3 == 0:
                    pfa.setFlux(float(i), 0.1*i)
                    pfa.setEllipticity(0.3*i, 0.4*i, 0.5*i,
                                       0.03*i, 0.04*i, 0.05*i)
                elif i % 3 == 1:
                    pfa.setFlux(float(i), None)
                    pfa.setEllipticity(0.3*i, 0.4*i, 0.5*i, None, None, None)
                else:
                    pfa.setFlux(None, None)
                    pfa.setEllipticity(None, None, None, None, None, None)
                sca.setPerFilterAttributes(pfa)
            self.clusters.append(sca)

    def tearDown(self):
        del self.clusters

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
        inp = cluster.PersistableSourceClusterVector(self.clusters)
        pol = policy.Policy()
        p = persistence.Persistence.getPersistence(pol)
        dp = base.PropertySet()
        psl = persistence.StorageList()
        rsl = persistence.StorageList()
        f = tempfile.NamedTemporaryFile()
        try:
            loc  = persistence.LogicalLocation(f.name)
            psl.append(p.getPersistStorage("BoostStorage", loc))
            p.persist(inp, psl, dp)
            rsl.append(p.getRetrieveStorage("BoostStorage", loc))
            persistable = p.unsafeRetrieve(
                "PersistableSourceClusterVector", rsl, dp)
            outp = cluster.PersistableSourceClusterVector.swigConvert(persistable)
            clusters = outp.getClusters()
            self.assertEqual(len(self.clusters), len(clusters))
            for sc1, sc2 in izip(self.clusters, clusters):
                self.assertEqual(sc1, sc2)
        finally:
            f.close()

    def testDatabasePersistence(self):
        """Tests database persistence of
        lsst::ap::cluster::SourceClusterAttributes
        """
        if not persistence.DbAuth.available("lsst10.ncsa.uiuc.edu", "3306"):
            # skip database test
            return

        # generate a (hopefully unique) table name for test data
        tableName = platform.node().replace('.', '_') + str(os.getpid())

        # register filters expected by the formatter
        image.Filter.reset()
        image.Filter.define(image.FilterProperty("u"), 0)
        image.Filter.define(image.FilterProperty("g"), 1)  
        image.Filter.define(image.FilterProperty("r"), 2)  
        image.Filter.define(image.FilterProperty("i"), 3)  
        image.Filter.define(image.FilterProperty("z"), 4)  
        image.Filter.define(image.FilterProperty("y"), 5)  

        # setup persistence machinery
        inp = cluster.PersistableSourceClusterVector(self.clusters)
        pol = policy.Policy()
        p = persistence.Persistence.getPersistence(pol)
        dp = base.PropertySet()
        pol.set("Formatter.PersistableSourceClusterVector.Clusters.templateTableName", "Object")
        pol.set("Formatter.PersistableSourceClusterVector.Clusters.tableNamePattern", tableName)
        loc = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test_object_pt1")
        dp.setString("itemName", "Clusters")
        for storageType in ("DbStorage", "DbTsvStorage"):
            psl = persistence.StorageList()
            rsl = persistence.StorageList()
            psl.append(p.getPersistStorage(storageType, loc))
            rsl.append(p.getRetrieveStorage(storageType, loc))
            try:
                p.persist(inp, psl, dp)
                persistable = p.unsafeRetrieve(
                    "PersistableSourceClusterVector", rsl, dp)
                outp = cluster.PersistableSourceClusterVector.swigConvert(persistable)
                clusters = outp.getClusters()
                self.assertEqual(len(self.clusters), len(clusters))
                for sc1, sc2 in izip(self.clusters, clusters):
                    # cannot test for exact equality since the database formatter
                    # does unit conversion
                    self.assertEqual(sc1.getClusterId(), sc2.getClusterId())
                    self.assertEqual(sc1.getNumObs(), sc2.getNumObs())
                    self.assertAlmostEqual(sc1.getRa(), sc2.getRa(), 14)
                    self.assertAlmostEqual(sc1.getDec(), sc2.getDec(), 14)
                    if sc1.getRaSigma() is None:
                        self.assertEqual(sc2.getRaSigma(), None)
                    else:
                        self.assertAlmostEqual(sc1.getRaSigma(), sc2.getRaSigma(), 6)
                    if sc1.getDecSigma() is None:
                        self.assertEqual(sc2.getDecSigma(), None)
                    else:
                        self.assertAlmostEqual(sc1.getDecSigma(), sc2.getDecSigma(), 6)
                    if sc1.getRaDecCov() is None:
                        self.assertEqual(sc2.getRaDecCov(), None)
                    else:
                        self.assertAlmostEqual(sc1.getRaDecCov(), sc2.getRaDecCov())
                    self.assertEqual(sc1.getPerFilterAttributes(),
                                     sc2.getPerFilterAttributes())
            finally:
                db = persistence.DbStorage()
                db.setPersistLocation(loc)
                db.startTransaction()
                try:
                    db.dropTable(tableName)
                except:
                    pass
                db.endTransaction()


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

