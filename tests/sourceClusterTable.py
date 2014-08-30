#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2012 LSST Corporation.
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

"""
Tests for lsst.ap.cluster.SourceClusterTable

Run with:
   python sourceClusterTable.py
or
   python
   >>> import sourceClusterTable; sourceClusterTable.run()
"""

import sys
import os
import unittest
import math
import numpy

import lsst.utils.tests
import lsst.pex.exceptions
import lsst.afw.table
import lsst.afw.geom
import lsst.afw.coord
import lsst.ap.cluster
import lsst.ap.utils


def makeCov(size, dtype):
    m = numpy.array(numpy.random.randn(size, size), dtype=dtype)
    return numpy.dot(m, m.transpose())

class SourceClusterTableTestCase(unittest.TestCase):

    def fillRecord(self, record):
        decl = numpy.random.randn()
        while decl <= -0.5*math.pi or decl >= 0.5*math.pi:
            decl = numpy.random.randn()
        record.setCoord(lsst.afw.coord.IcrsCoord(
            numpy.random.randn() * lsst.afw.geom.radians,
            decl * lsst.afw.geom.radians,
        ))
        record.set(self.coordErrKey, makeCov(2, numpy.float32))
        record.set(self.numSourcesKey, numpy.random.randint(0, 10000))
        record.set(self.timeMinKey, numpy.random.randn())
        record.set(self.timeMaxKey, numpy.random.randn())
        record.set(self.uFluxKey, numpy.random.randn())
        record.set(self.uFluxErrKey, numpy.random.randn())
        record.set(self.uFluxCountKey, numpy.random.randint(0, 10000))
        record.set(self.zFluxKey, numpy.random.randn())
        record.set(self.zFluxErrKey, numpy.random.randn())
        record.set(self.zFluxCountKey, numpy.random.randint(0, 10000))
        record.set(self.gShapeKey, lsst.afw.geom.ellipses.Quadrupole(*numpy.random.randn(3)))
        record.set(self.gShapeErrKey, makeCov(3, numpy.float32))
        record.set(self.gShapeCountKey, numpy.random.randint(0, 10000))
        record.set(self.rNumSourcesKey, numpy.random.randint(0, 10000))
        record.set(self.rTimeMinKey, numpy.random.randn())
        record.set(self.rTimeMaxKey, numpy.random.randn())
        record.set(self.flag1Key, numpy.random.randn() > 0)
        record.set(self.flag2Key, numpy.random.randn() > 0)

    def setUp(self):
        self.schema = lsst.ap.cluster.SourceClusterTable.makeMinimalSchema()
        self.schema.setVersion(0)
        self.coordErrKey = self.schema.addField("coord.err", type="CovPointF")
        self.numSourcesKey = self.schema.addField("obs.num", type="I")
        self.timeMinKey = self.schema.addField("obs.time.min", type="D")
        self.timeMaxKey = self.schema.addField("obs.time.max", type="D")
        self.uFluxKey = self.schema.addField("u.flux", type="D")
        self.uFluxErrKey = self.schema.addField("u.flux.err", type="D")
        self.uFluxCountKey = self.schema.addField("u.flux.count", type="I")
        self.zFluxKey = self.schema.addField("z.flux", type="D")
        self.zFluxErrKey = self.schema.addField("z.flux.err", type="D")
        self.zFluxCountKey = self.schema.addField("z.flux.count", type="I")
        self.gShapeKey = self.schema.addField("g.shape", type="MomentsD")
        self.gShapeErrKey = self.schema.addField("g.shape.err", type="CovMomentsF")
        self.gShapeCountKey = self.schema.addField("g.shape.count", type="I")
        self.rNumSourcesKey = self.schema.addField("r.obs.num", type="I")
        self.rTimeMinKey = self.schema.addField("r.obs.time.min", type="D")
        self.rTimeMaxKey = self.schema.addField("r.obs.time.max", type="D")
        self.flag1Key = self.schema.addField("flag1", type="Flag")
        self.flag2Key = self.schema.addField("flag2", type="Flag")
        self.idFactory = lsst.ap.cluster.SourceClusterIdFactory(1) 
        self.table = lsst.ap.cluster.SourceClusterTable.make(self.schema, self.idFactory)
        self.catalog = lsst.ap.cluster.SourceClusterCatalog(self.table)
        self.record = self.catalog.addNew()
        self.fillRecord(self.record)
        self.record.setId((1<<32) + 50)
        self.fillRecord(self.catalog.addNew())
        self.fillRecord(self.catalog.addNew())

    def tearDown(self):
        del self.schema
        del self.record
        del self.idFactory
        del self.table
        del self.catalog

    def checkCanonical(self):
        self.assertEqual(self.table.getPsfFluxDefinition("u"), "u.flux")
        self.assertEqual(self.record.get(self.uFluxKey), self.record.getPsfFlux("u"))
        self.assertEqual(self.record.get(self.uFluxErrKey), self.record.getPsfFluxErr("u"))
        self.assertEqual(self.record.get(self.uFluxCountKey), self.record.getPsfFluxCount("u"))
        self.assertEqual(self.table.getModelFluxDefinition("z"), "z.flux")
        self.assertEqual(self.record.get(self.zFluxKey), self.record.getModelFlux("z"))
        self.assertEqual(self.record.get(self.zFluxErrKey), self.record.getModelFluxErr("z"))
        self.assertEqual(self.record.get(self.zFluxCountKey), self.record.getModelFluxCount("z"))
        self.assertEqual(self.table.getShapeDefinition("g"), "g.shape")
        self.assertEqual(self.record.get(self.gShapeKey), self.record.getShape("g"))
        self.assert_(numpy.all(self.record.get(self.gShapeErrKey) == self.record.getShapeErr("g")))
        self.assertEqual(self.record.get(self.gShapeCountKey), self.record.getShapeCount("g"))

    def testCanonical1(self):
        self.table.defineCoordErr(self.coordErrKey)
        self.table.defineNumSources(self.numSourcesKey)
        self.table.defineTimeMin(self.timeMinKey)
        self.table.defineTimeMax(self.timeMaxKey)
        self.table.defineNumSources("r", self.rNumSourcesKey)
        self.table.defineTimeMin("r", self.rTimeMinKey)
        self.table.defineTimeMax("r", self.rTimeMaxKey)
        self.assertEqual(self.table.getCoordErrDefinition(), "coord.err")
        self.assert_(numpy.all(self.record.get(self.coordErrKey) == self.record.getCoordErr()))
        self.assertEqual(self.table.getNumSourcesDefinition(), "obs.num")
        self.assertEqual(self.record.get(self.numSourcesKey), self.record.getNumSources())
        self.assertEqual(self.table.getTimeMinDefinition(), "obs.time.min")
        self.assertEqual(self.record.get(self.timeMinKey), self.record.getTimeMin())
        self.assertEqual(self.table.getTimeMaxDefinition(), "obs.time.max")
        self.assertEqual(self.record.get(self.timeMaxKey), self.record.getTimeMax())
        self.assertEqual(self.table.getNumSourcesDefinition("r"), "r.obs.num")
        self.assertEqual(self.record.get(self.rNumSourcesKey), self.record.getNumSources("r"))
        self.assertEqual(self.table.getTimeMinDefinition("r"), "r.obs.time.min")
        self.assertEqual(self.record.get(self.rTimeMinKey), self.record.getTimeMin("r"))
        self.assertEqual(self.table.getTimeMaxDefinition("r"), "r.obs.time.max")
        self.assertEqual(self.record.get(self.rTimeMaxKey), self.record.getTimeMax("r"))
        self.table.definePsfFlux("u", self.uFluxKey, self.uFluxErrKey, self.uFluxCountKey)
        self.table.defineModelFlux("z", self.zFluxKey, self.zFluxErrKey, self.zFluxCountKey)
        self.table.defineShape("g", self.gShapeKey, self.gShapeErrKey, self.gShapeCountKey)
        self.checkCanonical()

    def testCanonical2(self):
        self.table.definePsfFlux("u", "flux")
        self.table.defineModelFlux("z", "flux")
        self.table.defineShape("g", "shape")
        self.checkCanonical()

    def testSorting(self):
        self.assertFalse(self.catalog.isSorted())
        self.catalog.sort()
        self.assert_(self.catalog.isSorted())
        r = self.catalog.find((1<<32) + 2)
        self.assertEqual(r["id"], (1<<32) + 2)
        r = self.catalog.find(500)
        self.assert_(r is None)

    def testConversion(self):
        catalog1 = self.catalog.cast(lsst.ap.cluster.SourceClusterCatalog)
        catalog2 = self.catalog.cast(lsst.afw.table.SimpleCatalog)
        catalog3 = self.catalog.cast(lsst.ap.cluster.SourceClusterCatalog, deep=True)
        catalog4 = self.catalog.cast(lsst.afw.table.SimpleCatalog, deep=True)
        self.assertEqual(self.catalog.table, catalog1.table)
        self.assertEqual(self.catalog.table, catalog2.table)
        self.assertNotEqual(self.catalog.table, catalog3.table)
        self.assertNotEqual(self.catalog.table, catalog3.table)
        for r, r1, r2, r3, r4 in zip(self.catalog, catalog1, catalog2, catalog3, catalog4):
            self.assertEqual(r, r1)
            self.assertEqual(r, r2)
            self.assertNotEqual(r, r3)
            self.assertNotEqual(r, r4)
            self.assertEqual(r.getId(), r3.getId())
            self.assertEqual(r.getId(), r4.getId())

    def testColumnView(self):
        cols1 = self.catalog.getColumnView()
        cols2 = self.catalog.columns
        self.assert_(cols1 is cols2)
        self.assert_(isinstance(cols1, lsst.ap.cluster.SourceClusterColumnView))
        self.table.definePsfFlux("u", "flux")
        self.table.defineModelFlux("z", "flux")
        self.table.defineShape("g", "shape")
        self.assert_((cols2["u.flux"] == cols2.getPsfFlux("u")).all())
        self.assert_((cols2["z.flux"] == cols2.getModelFlux("z")).all())
        self.assert_((cols2["g.shape.xx"] == cols2.getIxx("g")).all())
        self.assert_((cols2["g.shape.yy"] == cols2.getIyy("g")).all())
        self.assert_((cols2["g.shape.xy"] == cols2.getIxy("g")).all())

    def testForwarding(self):
        """Verify that Catalog forwards unknown methods to its table and/or columns."""
        self.table.definePsfFlux("u", "flux")
        self.table.defineModelFlux("z", "flux")
        self.table.defineShape("g", "shape")
        self.assert_((self.catalog.columns["u.flux"] == self.catalog["u.flux"]).all())
        self.assert_((self.catalog.columns[self.uFluxKey] == self.catalog.get(self.uFluxKey)).all())
        self.assert_((self.catalog.columns.get(self.uFluxKey) == self.catalog.getPsfFlux("u")).all())
        self.assertEqual(self.uFluxKey, self.catalog.getPsfFluxKey("u"))
        self.assertRaises(AttributeError, lambda c: c.foo(), self.catalog)

    def testBitsColumn(self):
        allBits = self.catalog.getBits()
        someBits = self.catalog.getBits(["flag2"])
        self.assertEqual(allBits.getMask("flag1"), 0x1)
        self.assertEqual(allBits.getMask("flag2"), 0x2)
        self.assertEqual(someBits.getMask(self.flag2Key), 0x1)
        self.assert_(((allBits.array & 0x1 != 0) == self.catalog.columns["flag1"]).all())
        self.assert_(((allBits.array & 0x2 != 0) == self.catalog.columns["flag2"]).all())
        self.assert_(((someBits.array & 0x1 != 0) == self.catalog.columns["flag2"]).all())

    def testCast(self):
        baseCat = self.catalog.cast(lsst.afw.table.BaseCatalog)
        clusterCat = baseCat.cast(lsst.ap.cluster.SourceClusterCatalog)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    lsst.utils.tests.init()
    suites = []
    suites += unittest.makeSuite(SourceClusterTableTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
