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
import lsst.ap.match as match


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
            sca.setObsTime(float(i), float(i + 1), i + 0.5)
            if i % 4 == 0:
                sca.setPsPosition(0.2*i, 0.1*(i - 10), 0.01*i, 0.05*i, 0.03*i)
                sca.setSgPosition(0.2*i, 0.1*(i - 10), 0.01*i, 0.05*i, 0.03*i)
            elif i % 4 == 1:
                sca.setPsPosition(0.2*i, 0.1*(i - 10), 0.01*i, 0.05*i, None)
                sca.setSgPosition(0.2*i, 0.1*(i - 10), 0.01*i, 0.05*i, None)
            elif i % 4 == 2:
                sca.setPsPosition(0.2*i, 0.1*(i - 10), None, None, None)
                sca.setSgPosition(0.2*i, 0.1*(i - 10), None, None, None)
            else:
                sca.setPsPosition(0.2*i, 0.1*(i - 10), None, None, None)
                sca.setSgPosition(None, None, None, None, None)
            # add in per-filter attributes
            for j in xrange(i % 7):
                pfa = cluster.PerFilterSourceClusterAttributes()
                pfa.setFilterId(j)
                pfa.setNumObs(j + 1)
                sca.setNumObs(sca.getNumObs() + pfa.getNumObs())
                pfa.setNumPsFluxSamples(max(0, j - 2))
                pfa.setNumSgFluxSamples(max(0, j - 2))
                pfa.setNumGaussianFluxSamples(max(0, j - 2))
                if i % 3 == 0:
                    pfa.setPsFlux(float(i), 0.1*i)
                    pfa.setSgFlux(float(i), 0.1*i)
                    pfa.setGaussianFlux(float(i), 0.1*i)
                    pfa.setEllipticity(0.3*i, 0.4*i, 0.5*i,
                                       0.03*i, 0.04*i, 0.05*i)
                elif i % 3 == 1:
                    pfa.setPsFlux(float(i), None)
                    pfa.setSgFlux(float(i), None)
                    pfa.setGaussianFlux(float(i), None)
                    pfa.setEllipticity(0.3*i, 0.4*i, 0.5*i, None, None, None)
                else:
                    pfa.setPsFlux(None, None)
                    pfa.setSgFlux(None, None)
                    pfa.setGaussianFlux(None, None)
                    pfa.setEllipticity(None, None, None, None, None, None)
                sca.setPerFilterAttributes(pfa)
            self.clusters.append(sca)
        self.exposures = match.ExposureInfoMap()
        for filterId, filter in enumerate("ugrizy"):
            ps = base.PropertySet()
            ps.setLong("scienceCcdExposureId", filterId)
            ps.setString("FILTER", filter)
            ps.setString("TIME-MID", "2000-01-01T11:59:28.000000000Z")
            ps.setDouble("EXPTIME", 10.0)
            ps.setString("RADESYS", "FK5")
            ps.setDouble("EQUINOX", 2000.0)
            ps.setString("CTYPE1", "RA---TAN")
            ps.setString("CTYPE2", "DEC--TAN")
            ps.setString("CUNIT1", "deg")
            ps.setString("CUNIT2", "deg")
            ps.setInt("NAXIS1", 3600)
            ps.setInt("NAXIS2", 3600)
            ps.setDouble("CRPIX1", 1800.5)
            ps.setDouble("CRPIX2", 1800.5)
            ps.setDouble("CRVAL1", 0.0)
            ps.setDouble("CRVAL2", 0.0)
            ps.setDouble("CD1_1", 1.0/3600.0)
            ps.setDouble("CD1_2", 0.0)
            ps.setDouble("CD2_1", 0.0)
            ps.setDouble("CD2_2", 1.0/3600.0)
            ps.setDouble("FLUXMAG0", 1.0)
            ps.setDouble("FLUXMAG0ERR", 0.0)
            self.exposures.insert(match.ExposureInfo(ps))

    def tearDown(self):
        del self.clusters
        del self.exposures

    def testPosition(self):
        """Tests cluster position computation.
        """
        sources = detection.SourceSet()
        for i in xrange(5):
            s = detection.Source()
            s.setAmpExposureId(0)
            s.setFilterId(0)
            s.setRa(0.1*i)
            s.setDec(0.0)
            sources.append(s)
            if i != 2:
                s = detection.Source()
                s.setRa(0.2)
                s.setDec(0.1*(i - 2))
                sources.append(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 0, 0, 0, 0, 0)
        self.assertAlmostEqual(sca.getRaPs(), 0.2)
        self.assertAlmostEqual(sca.getDecPs(), 0.0)
        self.assertAlmostEqual(sca.getRaPsSigma(), math.sqrt(0.1/72))
        self.assertAlmostEqual(sca.getDecPsSigma(), math.sqrt(0.1/72))
        self.assertAlmostEqual(sca.getRaDecPsCov(), 0.0)

    def testTimes(self):
        """Tests observation time range computation.
        """
        sources = detection.SourceSet()
        for i in xrange(6):
            s = detection.Source()
            s.setTaiMidPoint(i)
            s.setFilterId(i // 3)
            s.setAmpExposureId(i // 3)
            sources.append(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 0, 0, 0, 0, 0)
        self.assertEqual(sca.getNumObs(), 6)
        self.assertEqual(sca.getEarliestObsTime(), 0)
        self.assertEqual(sca.getLatestObsTime(), 5)
        self.assertEqual(sca.getMeanObsTime(), 2.5)
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
            cluster.SourceClusterAttributes.setObsTime,
            sca, 1.0, -1.0, 0.0)
        for m in (0.0, 3.0):
            self.assertRaises(
                exceptions.LsstException,
                cluster.SourceClusterAttributes.setObsTime,
                sca, 1.0, 2.0, m)

    def testEllipticity(self):
        """Tests computation of ellipticity parameters.
        """
        sources = detection.SourceSet()
        esources = []
        for i in xrange(1, 5):
            s = detection.Source()
            s.setFilterId(0)
            s.setAmpExposureId(0)
            s.setIxx(i)
            s.setIxxErr(0.1)
            s.setIyy(2.0*i)
            s.setIyyErr(0.1)
            s.setIxy(i*0.5)
            s.setIxyErr(0.1)
            if i == 4:
                s.setFlagForDetection(1)
            else:
                esources.append(_eparams(s))
            sources.append(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 0, 0, 0, 0, 1)
        pfa = sca.getPerFilterAttributes(0)
        self.assertEqual(pfa.getNumObs(), 4)
        self.assertEqual(pfa.getNumEllipticitySamples(), 3)
        eparams = []
        for i in xrange(3):
            eparams.append(sum(e[i] for e in esources) / len(esources))
        self.assertAlmostEqual(pfa.getE1(), eparams[0], 1)
        self.assertAlmostEqual(pfa.getE2(), eparams[1], 1)
        self.assertAlmostEqual(pfa.getRadius(), math.radians(eparams[2])/3600.0, 1)
        # Test with a single source
        sources.clear()
        s = detection.Source()
        s.setAmpExposureId(0)
        s.setFilterId(0)
        s.setIxx(1.0)
        s.setIxxErr(0.1)
        s.setIyy(2.0)
        s.setIyyErr(0.1)
        s.setIxy(1.0)
        s.setIxyErr(0.1)
        sources.append(s)
        eparams = _eparams(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 0, 0, 0, 0, 0)
        pfa = sca.getPerFilterAttributes(0)
        self.assertAlmostEqual(pfa.getE1(), eparams[0], 6)
        self.assertAlmostEqual(pfa.getE2(), eparams[1], 6)
        self.assertAlmostEqual(pfa.getRadius(), math.radians(eparams[2])/3600.0, 6)
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
            s.setFilterId(0)
            s.setAmpExposureId(0)
            s.setPsfFlux(i)
            s.setPsfFluxErr(1.0)
            s.setModelFlux(i)
            s.setModelFluxErr(2.0)
            s.setInstFlux(i)
            s.setInstFluxErr(4.0)
            s.setFlagForDetection(int(i == 5))
            sources.append(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 1, 1, 0, 1, 0)
        pfa = sca.getPerFilterAttributes(0)
        self.assertEqual(pfa.getNumObs(), 6)
        self.assertEqual(pfa.getNumPsFluxSamples(), 5)
        self.assertEqual(pfa.getPsFlux(), 2)
        self.assertAlmostEqual(pfa.getPsFluxSigma(), 0.5 * math.sqrt(2.0), 6)
        self.assertEqual(pfa.getSgFlux(), 2)
        self.assertAlmostEqual(pfa.getSgFluxSigma(), 0.5 * math.sqrt(2.0), 6)
        self.assertEqual(pfa.getGaussianFlux(), 2)
        self.assertAlmostEqual(pfa.getGaussianFluxSigma(), 0.5 * math.sqrt(2.0), 6)

        # Test with a single source
        sources.clear()
        s = detection.Source()
        s.setFilterId(0)
        s.setAmpExposureId(0)
        s.setPsfFlux(1.0)
        s.setPsfFluxErr(1.0)
        s.setModelFlux(2.0)
        s.setModelFluxErr(2.0)
        s.setInstFlux(4.0)
        s.setInstFluxErr(4.0)
        sources.append(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 1, 1, 0, 1, 0)
        pfa = sca.getPerFilterAttributes(0)
        self.assertEqual(pfa.getPsFlux(), 1.0)
        self.assertEqual(pfa.getPsFluxSigma(), 1.0)
        self.assertEqual(pfa.getNumPsFluxSamples(), 1)
        self.assertEqual(pfa.getSgFlux(), 2.0)
        self.assertEqual(pfa.getSgFluxSigma(), 2.0)
        self.assertEqual(pfa.getNumSgFluxSamples(), 1)
        self.assertEqual(pfa.getGaussianFlux(), 4.0)
        self.assertEqual(pfa.getGaussianFluxSigma(), 4.0)
        self.assertEqual(pfa.getNumGaussianFluxSamples(), 1)
        # Test with bad fluxes
        sources.clear()
        for i in xrange(6):
            s = detection.Source()
            s.setFilterId(0)
            s.setAmpExposureId(0)
            s.setPsfFlux(i)
            s.setFlagForDetection(1)
            sources.append(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 1, 1, 0, 1, 0)
        pfa = sca.getPerFilterAttributes(0)
        self.assertEqual(pfa.getNumPsFluxSamples(), 0)
        self.assertEqual(pfa.getNumSgFluxSamples(), 0)
        self.assertEqual(pfa.getNumGaussianFluxSamples(), 0)
        self.assertEqual(pfa.getPsFlux(), None)
        self.assertEqual(pfa.getPsFluxSigma(), None)
        self.assertEqual(pfa.getSgFlux(), None)
        self.assertEqual(pfa.getSgFluxSigma(), None)
        self.assertEqual(pfa.getGaussianFlux(), None)
        self.assertEqual(pfa.getGaussianFluxSigma(), None)
        self.assertEqual(pfa.getNumObs(), 6)
        # Test error checking
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setPsFlux,
                          pfa, None, 1.0)
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setPsFlux,
                          pfa, 1.0, -1.0)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumPsFluxSamples,
            pfa, -1)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumPsFluxSamples,
            pfa, cluster.PerFilterSourceClusterAttributes.NSAMPLE_MASK + 1)
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setSgFlux,
                          pfa, None, 1.0)
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setSgFlux,
                          pfa, 1.0, -1.0)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumSgFluxSamples,
            pfa, -1)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumSgFluxSamples,
            pfa, cluster.PerFilterSourceClusterAttributes.NSAMPLE_MASK + 1)
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setGaussianFlux,
                          pfa, None, 1.0)
        self.assertRaises(exceptions.LsstException,
                          cluster.PerFilterSourceClusterAttributes.setGaussianFlux,
                          pfa, 1.0, -1.0)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumGaussianFluxSamples,
            pfa, -1)
        self.assertRaises(
            exceptions.LsstException,
            cluster.PerFilterSourceClusterAttributes.setNumGaussianFluxSamples,
            pfa, cluster.PerFilterSourceClusterAttributes.NSAMPLE_MASK + 1)

    def testSourceClusterAttributes(self):
        """Tests features of the Python wrapper for
        lsst::ap::cluster::SourceClusterAttributes
        """
        sources = detection.SourceSet()
        eparams = []
        for i in xrange(6):
            s = detection.Source()
            s.setAmpExposureId(i)
            s.setFilterId(i)
            s.setPsfFlux(i)
            s.setPsfFluxErr(i + 1)
            s.setModelFlux(2*i)
            s.setModelFluxErr(2*i + 1)
            s.setInstFlux(3*i)
            s.setInstFluxErr(3*i + 1)
            s.setIxx(i + 1)
            s.setIxxErr(0.1)
            s.setIyy(2.0*i + 2.0)
            s.setIyyErr(0.1)
            s.setIxy(0.5*i + 0.5)
            s.setIxyErr(0.1)
            eparams.append(_eparams(s))
            sources.append(s)
        sca = cluster.SourceClusterAttributes(0)
        sca.computeAttributes(sources, self.exposures, 1.0, 0, 0, 0, 0, 0)
        filters = sorted(sca.getFilterIds())
        self.assertEqual(filters, range(6))
        filterMap = sca.getPerFilterAttributes()
        for f in filters:
            self.assertTrue(sca.hasFilter(f))
            self.assertTrue(f in filterMap)
            pfa = sca.getPerFilterAttributes(f)
            self.assertEqual(pfa.getNumObs(), 1)
            self.assertEqual(pfa.getNumPsFluxSamples(), 1)
            self.assertEqual(pfa.getNumSgFluxSamples(), 1)
            self.assertEqual(pfa.getNumGaussianFluxSamples(), 1)
            self.assertEqual(pfa.getNumEllipticitySamples(), 1)
            self.assertEqual(pfa.getPsFlux(), f)
            self.assertEqual(pfa.getPsFluxSigma(), f + 1)
            self.assertEqual(pfa.getSgFlux(), 2*f)
            self.assertEqual(pfa.getSgFluxSigma(), 2*f + 1)
            self.assertEqual(pfa.getGaussianFlux(), 3*f)
            self.assertEqual(pfa.getGaussianFluxSigma(), 3*f + 1)
            self.assertAlmostEqual(pfa.getE1(), eparams[f][0], 6)
            self.assertAlmostEqual(pfa.getE2(), eparams[f][1], 6)
            self.assertAlmostEqual(pfa.getRadius(),
                                   math.radians(eparams[f][2])/3600.0, 6)
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
        for f, fid in zip("ugrizy", range(6)):
            image.Filter.define(image.FilterProperty(f), fid)

        # setup persistence machinery
        inp = cluster.PersistableSourceClusterVector(self.clusters)
        pol = policy.Policy()
        p = persistence.Persistence.getPersistence(pol)
        dp = base.PropertySet()
        pol.set("Formatter.PersistableSourceClusterVector.Clusters.templateTableName", "Object")
        pol.set("Formatter.PersistableSourceClusterVector.Clusters.tableNamePattern", tableName)
        loc = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test_object_pt1_2")
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
                    self.assertAlmostEqual(sc1.getRaPs(), sc2.getRaPs(), 14)
                    self.assertAlmostEqual(sc1.getDecPs(), sc2.getDecPs(), 14)
                    if sc1.getRaPsSigma() == None:
                        self.assertEqual(sc2.getRaPsSigma(), None)
                    else:
                        self.assertAlmostEqual(sc1.getRaPsSigma(), sc2.getRaPsSigma(), 6)
                    if sc1.getDecPsSigma() == None:
                        self.assertEqual(sc2.getDecPsSigma(), None)
                    else:
                        self.assertAlmostEqual(sc1.getDecPsSigma(), sc2.getDecPsSigma(), 6)
                    if sc1.getRaDecPsCov() == None:
                        self.assertEqual(sc2.getRaDecPsCov(), None)
                    else:
                        self.assertAlmostEqual(sc1.getRaDecPsCov(), sc2.getRaDecPsCov(), 6)
                    if sc1.getRaSg() == None:
                        self.assertEqual(sc2.getRaSg(), None)
                    else:
                        self.assertAlmostEqual(sc1.getRaSg(), sc2.getRaSg(), 14)
                    if sc1.getDecSg() == None:
                        self.assertEqual(sc2.getDecSg(), None)
                    else:
                        self.assertAlmostEqual(sc1.getDecSg(), sc2.getDecSg(), 14)
                    if sc1.getRaSgSigma() == None:
                        self.assertEqual(sc2.getRaSgSigma(), None)
                    else:
                        self.assertAlmostEqual(sc1.getRaSgSigma(), sc2.getRaSgSigma(), 6)
                    if sc1.getDecSgSigma() == None:
                        self.assertEqual(sc2.getDecSgSigma(), None)
                    else:
                        self.assertAlmostEqual(sc1.getDecSgSigma(), sc2.getDecSgSigma(), 6)
                    if sc1.getRaDecSgCov() == None:
                        self.assertEqual(sc2.getRaDecSgCov(), None)
                    else:
                        self.assertAlmostEqual(sc1.getRaDecSgCov(), sc2.getRaDecSgCov(), 6)
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

