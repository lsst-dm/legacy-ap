# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2012 LSST Corporation.
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
import numpy
import os
import pdb
import unittest

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.image.utils as afwImageUtils
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwGeomEllipses
import lsst.afw.table as afwTable
import lsst.ap.cluster as apCluster
import lsst.ap.match as apMatch


# Usually defined by lsst.obs.lsstSim.LsstSimMapper
afwImageUtils.defineFilter('u', 364.59)
afwImageUtils.defineFilter('g', 476.31)
afwImageUtils.defineFilter('r', 619.42)
afwImageUtils.defineFilter('i', 752.06)
afwImageUtils.defineFilter('z', 866.85)
afwImageUtils.defineFilter('y', 971.68)

class SourceClusterAttributesTestCase(unittest.TestCase):
    """Test case for source cluster attribute computation.
    """
    def setUp(self):
        # Setup bogus exposure information
        self.exposures = apMatch.ExposureInfoMap()
        for filter in "ugrizy":
            filterId = afwImage.Filter(filter).getId()
            ps = dafBase.PropertySet()
            ps.setLong("Computed_ccdExposureId", filterId)
            ps.setString("FILTER", filter)
            ps.setString("TIME-MID", "2000-01-%02dT11:59:28.000000000Z" % (filterId + 1))
            ps.setDouble("EXPTIME", 10.0)
            ps.setString("RADESYS", "ICRS")
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
            ps.setDouble("FLUXMAG0", 1.0)    # make calibration a no-op
            ps.setDouble("FLUXMAG0ERR", 0.0) # no error from calibration
            self.exposures.insert(apMatch.ExposureInfo(ps))
        # setup Source table/schema
        schema = afwTable.SourceTable.makeMinimalSchema()
        schema.addField("exposure.id", type="I8")
        schema.addField("exposure.filter.id", type="I4")
        schema.addField("exposure.time", type="F8")
        schema.addField("exposure.time.mid", type="F8")
        schema.addField("cluster.id", type="I8")
        schema.addField("cluster.coord", type="Coord")
        schema.addField("centroid", type="Point<F8>")
        schema.addField("centroid.err", type="Cov<Point<F8>>")
        schema.addField("centroid.flags", type="Flag")
        schema.addField("flux", type="F8")
        schema.addField("flux.err", type="F8")
        schema.addField("flux.flags", type="Flag")
        schema.addField("shape", type="Moments<F8>")
        schema.addField("shape.err", type="Cov<Moments<F8>>")
        schema.addField("shape.flags", type="Flag")
        schema.addField("flags.badflux", type="Flag")
        schema.addField("flags.badshape", type="Flag")
        self.sourceTable = afwTable.SourceTable.make(schema)
        self.sourceTable.defineCentroid("centroid")
        self.sourceTable.definePsfFlux("flux")
        self.sourceTable.defineShape("shape")
        # setup SourceCluster table/schemaa
        self.spConfig = apCluster.SourceProcessingConfig()
        self.spConfig.fluxFields = ["flux"]
        self.spConfig.shapeFields = ["shape"]
        self.clusterTable = apCluster.makeSourceClusterTable(
            self.sourceTable,
            apCluster.SourceClusterIdFactory(0),
            self.spConfig.makeControl())

    def tearDown(self):
        del self.clusterTable
        del self.spConfig
        del self.sourceTable
        del self.exposures

    def testBasic(self):
        """Test basic cluster attribute computation.
        """
        sources = afwTable.SourceCatalog(self.sourceTable)
        clusters = apCluster.SourceClusterCatalog(self.clusterTable)
        centroidKey = self.sourceTable.getCentroidKey()
        centroidErrKey = self.sourceTable.getCentroidErrKey()
        for i in xrange(5):
            s = sources.addNew()
            s["exposure.id"] = i
            s["exposure.filter.id"] = i
            s["exposure.time.mid"] = self.exposures.get(i).getEpoch()
            s.setRa((0.1 * (i - 2)) * afwGeom.radians)
            s.setDec(0.0 * afwGeom.radians)
            # note: test centroids are all for (ra,dec) = (0,0), not
            # the coordinates above. That's OK, so long as the
            # unweighted mean coordinates are (0,0)
            s.set(centroidKey, afwGeom.Point2D(1799.5, 1799.5))
            s.set(centroidErrKey, numpy.matrix([[10.0, 0.0], [0.0, 20.0]]))
            s.setRa((0.1 * (i - 2)) * afwGeom.radians)
            s.setDec(0.0 * afwGeom.radians)
            if i != 2:
                s = sources.addNew()
                s["exposure.id"] = i
                s["exposure.filter.id"] = i
                s["exposure.time.mid"] = self.exposures.get(i).getEpoch()
                s.set(centroidKey, afwGeom.Point2D(1799.5, 1799.5))
                s.set(centroidErrKey, numpy.matrix([[10.0, 0.0], [0.0, 20.0]]))
                s.setRa(0.0 * afwGeom.radians)
                s.setDec(0.1*(i - 2) * afwGeom.radians)
        cluster = clusters.addNew()
        seVec = apCluster.computeBasicAttributes(
            cluster, sources, self.exposures, self.spConfig.exposurePrefix)
        # make sure each source was associated with the correct exposure
        self.assertEqual(len(seVec), 9)
        for se in seVec:
            self.assertEqual(se.getSource()["exposure.id"], se.getExposureInfo().getId())
        for i in xrange(len(seVec) - 1):
            self.assertTrue(seVec[i].getSource()["exposure.filter.id"] <=
                            seVec[i + 1].getSource()["exposure.filter.id"])
        # test whether unweighted mean position was computed as expected
        coord = cluster.getCoord()
        self.assertAlmostEqual(coord.getLongitude().asRadians(), 0.0)
        self.assertAlmostEqual(coord.getLatitude().asRadians(), 0.0)
        coordErr = cluster.getCoordErr()
        self.assertAlmostEqual(coordErr[0,0], 0.1/8)
        self.assertAlmostEqual(coordErr[1,1], 0.1/8)
        self.assertEqual(coordErr[0,1], coordErr[1,0])
        self.assertAlmostEqual(coordErr[0,1], 0.0)
        # test whether weighted mean position was computed as expected
        self.assertEqual(cluster.getWeightedMeanCoordCount(), 9)
        coord = cluster.getWeightedMeanCoord()
        self.assertAlmostEqual(coord.getLongitude().asRadians(), 0.0)
        self.assertAlmostEqual(coord.getLatitude().asRadians(), 0.0)
        coordErr = cluster.getWeightedMeanCoordErr()
        self.assertEqual(coordErr[0,1], coordErr[1,0])
        self.assertAlmostEqual(coordErr[0,1], 0.0)
        scale2 = math.radians(1.0/3600.0)**2
        self.assertAlmostEqual(coordErr[0,0], scale2*10.0/9.0)
        self.assertAlmostEqual(coordErr[1,1], scale2*20.0/9.0)
        # test remaining basic attributes
        self.assertEqual(cluster["obs.count"], 9)
        self.assertEqual(cluster["obs.time.min"], self.exposures.get(0).getEpoch())
        self.assertEqual(cluster["obs.time.max"], self.exposures.get(4).getEpoch())
        self.assertAlmostEqual(cluster["obs.time.mean"], self.exposures.get(2).getEpoch())
        for i, filter in enumerate("ugriz"):
            self.assertEqual(cluster.getTimeMin(filter), cluster.getTimeMax(filter))
            self.assertEqual(cluster.getTimeMin(filter), self.exposures.get(i).getEpoch())
            self.assertEqual(cluster.getNumSources(filter), 1 if i == 2 else 2)


    def testFlux(self):
        """Test flux mean and error computation.
        """
        sources = afwTable.SourceCatalog(self.sourceTable)
        clusters = apCluster.SourceClusterCatalog(self.clusterTable)
        fluxKey = self.sourceTable.getPsfFluxKey()
        fluxErrKey = self.sourceTable.getPsfFluxErrKey()
        badfluxKey = self.sourceTable.getSchema().find("flags.badflux").key
        for i in xrange(6):
            s = sources.addNew()
            s["exposure.id"] = 0
            s["exposure.filter.id"] = 0
            s.setRa(0.0 * afwGeom.radians)
            s.setDec(0.0 * afwGeom.radians)
            s.set(fluxKey, i)
            s.set(fluxErrKey, 2.0)
            s.set(badfluxKey, (i == 5))
        cluster = clusters.addNew()
        seVec = apCluster.computeBasicAttributes(
            cluster, sources, self.exposures, self.spConfig.exposurePrefix)
        fkv = afwTable.FlagKeyVector()
        fkv.push_back(badfluxKey)
        apCluster.computeFluxMean(cluster, seVec, "flux", fkv, 1.0)
        self.assertEqual(cluster.getNumSources(), 6)
        self.assertEqual(cluster.getNumSources("u"), 6)
        self.assertEqual(cluster.getPsfFluxCount("u"), 5)
        self.assertAlmostEqual(cluster.getPsfFlux("u"), 2.0)
        self.assertAlmostEqual(cluster.getPsfFluxErr("u"), 0.5 * math.sqrt(2.0), 6)

        # Test with a single source
        s = sources[0]
        sources = afwTable.SourceCatalog(self.sourceTable)
        sources.append(s)
        cluster = clusters.addNew()
        seVec = apCluster.computeBasicAttributes(
            cluster, sources, self.exposures, self.spConfig.exposurePrefix)
        apCluster.computeFluxMean(cluster, seVec, "flux", fkv, 1.0)
        self.assertEqual(cluster.getNumSources(), 1)
        self.assertEqual(cluster.getNumSources("u"), 1)
        self.assertEqual(cluster.getPsfFluxCount("u"), 1)
        self.assertEqual(cluster.getPsfFlux("u"), 0.0)
        self.assertEqual(cluster.getPsfFluxErr("u"), 2.0)

        # Test with bad values
        sources = afwTable.SourceCatalog(self.sourceTable)
        for i in xrange(6):
            s = sources.addNew()
            s["exposure.id"] = 0
            s["exposure.filter.id"] = 0
            s.setRa(0.0 * afwGeom.radians)
            s.setDec(0.0 * afwGeom.radians)
            s.set(fluxKey, i)
            s.set(fluxErrKey, float('nan') if i % 1 == 0 else -1.0)
        cluster = clusters.addNew()
        seVec = apCluster.computeBasicAttributes(
            cluster, sources, self.exposures, self.spConfig.exposurePrefix)
        apCluster.computeFluxMean(cluster, seVec, "flux", fkv, 1.0)
        self.assertEqual(cluster.getNumSources(), 6)
        self.assertEqual(cluster.getNumSources("u"), 6)
        self.assertEqual(cluster.getPsfFluxCount("u"), 0)
        self.assertTrue(cluster.getPsfFlux("u") != cluster.getPsfFlux("u"))
        self.assertTrue(cluster.getPsfFluxErr("u") != cluster.getPsfFluxErr("u"))
        sources = afwTable.SourceCatalog(self.sourceTable)
        for i in xrange(6):
            s = sources.addNew()
            s["exposure.id"] = 0
            s["exposure.filter.id"] = 0
            s.setRa(0.0 * afwGeom.radians)
            s.setDec(0.0 * afwGeom.radians)
            s.set(fluxKey, float('nan') if i % 1 == 0 else 1.0)
            s.set(fluxErrKey, 1.0)
            s.set(badfluxKey, i % 1 != 0)
        cluster = clusters.addNew()
        seVec = apCluster.computeBasicAttributes(
            cluster, sources, self.exposures, self.spConfig.exposurePrefix)
        apCluster.computeFluxMean(cluster, seVec, "flux", fkv, 1.0)
        self.assertEqual(cluster.getNumSources(), 6)
        self.assertEqual(cluster.getNumSources("u"), 6)
        self.assertEqual(cluster.getPsfFluxCount("u"), 0)
        self.assertTrue(cluster.getPsfFlux("u") != cluster.getPsfFlux("u"))
        self.assertTrue(cluster.getPsfFluxErr("u") != cluster.getPsfFluxErr("u"))

    def testShape(self):
        """Tests computation of shape means.
        """
        sources = afwTable.SourceCatalog(self.sourceTable)
        clusters = apCluster.SourceClusterCatalog(self.clusterTable)
        centroidKey = self.sourceTable.getCentroidKey()
        shapeKey = self.sourceTable.getShapeKey()
        shapeErrKey = self.sourceTable.getShapeErrKey()
        badshapeKey = self.sourceTable.getSchema().find("flags.badshape").key
        for i in xrange(6):
            s = sources.addNew()
            s["exposure.id"] = 0
            s["exposure.filter.id"] = 0
            s.setRa(0.0 * afwGeom.radians)
            s.setDec(0.0 * afwGeom.radians)
            s.set(centroidKey, afwGeom.Point2D(1799.5, 1799.5))
            s.set(shapeKey, afwGeomEllipses.Quadrupole(i, 2.0*i, 0.5*i))
            s.set(shapeErrKey, numpy.identity(3))
            if i < 1 or i > 3:
                s.set(badshapeKey, True)
        cluster = clusters.addNew()
        seVec = apCluster.computeBasicAttributes(
            cluster, sources, self.exposures, self.spConfig.exposurePrefix)
        fkv = afwTable.FlagKeyVector()
        fkv.push_back(badshapeKey)
        apCluster.computeShapeMean(cluster, seVec, "shape", fkv)
        scale2 = math.radians(1.0/3600.0)**2
        self.assertEqual(cluster.getNumSources(), 6)
        self.assertEqual(cluster.getNumSources("u"), 6)
        self.assertEqual(cluster.getShapeCount("u"), 3)
        shape = cluster.getShape("u")
        self.assertAlmostEqual(shape.getIxx(), scale2 * 2.0)
        self.assertAlmostEqual(shape.getIxx(), scale2 * 4.0)
        self.assertAlmostEqual(shape.getIxx(), scale2 * 1.0)
        shapeErr = cluster.getShapeErr("u")
        self.assertAlmostEqual(shapeErr[0,1], 0.0)
        self.assertAlmostEqual(shapeErr[0,2], 0.0)
        self.assertAlmostEqual(shapeErr[1,2], 0.0)
        self.assertEqual(shapeErr[0,1], shapeErr[1,0])
        self.assertEqual(shapeErr[0,2], shapeErr[2,0])
        self.assertEqual(shapeErr[1,2], shapeErr[2,1])
        self.assertAlmostEqual(shapeErr[0,0], scale2**2/3)
        self.assertAlmostEqual(shapeErr[1,1], scale2**2/3)
        self.assertAlmostEqual(shapeErr[2,2], scale2**2/3)

        # Test with a single source
        s = sources[2]
        sources = afwTable.SourceCatalog(self.sourceTable)
        sources.append(s)
        cluster = clusters.addNew()
        seVec = apCluster.computeBasicAttributes(
            cluster, sources, self.exposures, self.spConfig.exposurePrefix)
        apCluster.computeShapeMean(cluster, seVec, "shape", fkv)
        self.assertEqual(cluster.getNumSources(), 1)
        self.assertEqual(cluster.getNumSources("u"), 1)
        self.assertEqual(cluster.getShapeCount("u"), 1)
        shape = cluster.getShape("u")
        self.assertAlmostEqual(shape.getIxx(), scale2 * 2.0)
        self.assertAlmostEqual(shape.getIxx(), scale2 * 4.0)
        self.assertAlmostEqual(shape.getIxx(), scale2 * 1.0)
        shapeErr = cluster.getShapeErr("u")
        self.assertAlmostEqual(shapeErr[0,1], 0.0)
        self.assertAlmostEqual(shapeErr[0,2], 0.0)
        self.assertAlmostEqual(shapeErr[1,2], 0.0)
        self.assertEqual(shapeErr[0,1], shapeErr[1,0])
        self.assertEqual(shapeErr[0,2], shapeErr[2,0])
        self.assertEqual(shapeErr[1,2], shapeErr[2,1])
        self.assertAlmostEqual(shapeErr[0,0], scale2**2)
        self.assertAlmostEqual(shapeErr[1,1], scale2**2)
        self.assertAlmostEqual(shapeErr[2,2], scale2**2)

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

