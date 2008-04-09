#! /usr/bin/env python

"""
Tests for C++ association pipeline result vector Python wrappers (including persistence)

Run with:
   python ResultTest.py
or
   python
   >>> import unittest; T=load("ResultTest"); unittest.TextTestRunner(verbosity=1).run(T.suite())
"""

import pdb
import unittest
import random
import time
import lsst.daf.base as dafBase
import lsst.pex.policy as policy
import lsst.daf.persistence as persistence
import lsst.utils.tests as tests
import lsst.ap.interface as ap


# ----------------------------------------------------------------
class MatchPairVecTestCase(unittest.TestCase):
    """A test case for MatchPair and MatchPairVec"""

    def setUp(self):
        n = 16 + random.randint(0,16)
        self.mpv1 = ap.MatchPairVec(n)
        self.mpv2 = ap.MatchPairVec()

        for i in xrange(n):
            self.mpv1[i].setFirst(i)
            mp = ap.MatchPair()
            mp.setFirst(i)
            mp.setSecond(i*20)
            mp.setDistance(random.random())
            self.mpv2.push_back(mp)

    def tearDown(self):
        del self.mpv1
        del self.mpv2

    def testIterable(self):
        """Tests that MatchPairVec instances can be iterated over"""
        j = 0
        for i in self.mpv1:
            assert i.getFirst() == j
            j += 1

    def testCopyAndCompare(self):
        """Tests copying, assignment and comparison of MatchPairVec instances"""
        mpv1Copy = ap.MatchPairVec(self.mpv1)
        mpv2Copy = ap.MatchPairVec(self.mpv2)
        assert mpv1Copy == self.mpv1
        assert mpv2Copy == self.mpv2
        mpv1Copy.swap(mpv2Copy)
        assert mpv1Copy == self.mpv2
        assert mpv2Copy == self.mpv1
        mpv1Copy.swap(mpv2Copy)
        if mpv1Copy.size() == 0:
            mpv1Copy.push_back(ap.MatchPair())
        else:
            mpv1Copy.pop_back()
        mp = ap.MatchPair()
        mp.setFirst(123456789876543210)
        mp.setSecond(876543210123456789)
        mp.setDistance(6.66)
        mpv2Copy.push_back(mp)
        assert mpv1Copy != self.mpv1
        assert mpv2Copy != self.mpv2

    def testInsertErase(self):
        """Tests inserting and erasing MatchPairVec elements"""
        mpv1Copy = ap.MatchPairVec(self.mpv1)
        mpv1Copy.insert(mpv1Copy.begin() + 5, ap.MatchPair())
        mpv1Copy.insert(mpv1Copy.begin() + 6, 4, ap.MatchPair())
        mpv1Copy.erase(mpv1Copy.begin() + 5)
        mpv1Copy.erase(mpv1Copy.begin() + 5, mpv1Copy.begin() + 9)
        assert mpv1Copy == self.mpv1

    def testSlice(self):
        slice = self.mpv1[1:5]
        j = 1
        for i in slice:
            print i
            assert i.getFirst() == j
            j += 1

    def testPersistence(self):
        if persistence.DbAuth.available():
            pol  = policy.PolicyPtr()
            pers = persistence.Persistence.getPersistence(pol)
            loc  = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test")
            dp   = dafBase.DataProperty.createPropertyNode("root")
            dp.addProperty(dafBase.DataProperty("visitId", int(time.clock())*16384 + random.randint(0,16383)))
            dp.addProperty(dafBase.DataProperty("itemName", "MatchPair"))
            stl = persistence.StorageList()
            stl.push_back(pers.getPersistStorage("DbStorage", loc))
            pers.persist(self.mpv1, stl, dp)
            stl = persistence.StorageList()
            stl.push_back(pers.getRetrieveStorage("DbStorage", loc))
            res = ap.MatchPairVecSharedPtr(pers.retrieve("MatchPairVector", stl, dp))
            db = persistence.DbStorage()
            db.setPersistLocation(loc)
            db.startTransaction()
            db.dropTable(ap.getTableName(pol, dp))
            db.endTransaction()
            assert(res == self.mpv1)


# ----------------------------------------------------------------
class IdPairVecTestCase(unittest.TestCase):
    """A test case for IdPair and IdPairVec"""

    def setUp(self):
        n = 16 + random.randint(0,16)
        self.ipv1 = ap.IdPairVec(n)
        self.ipv2 = ap.IdPairVec()

        for i in xrange(n):
            self.ipv1[i] = (i, 35)
            self.ipv2.push_back((i, -i))

    def tearDown(self):
        del self.ipv1
        del self.ipv2

    def testIterable(self):
        """Tests that IdPairVec instances can be iterated over"""
        j = 0
        for i in self.ipv1:
            assert i[0] == j
            j += 1

    def testCopyAndCompare(self):
        """Tests copying, assignment and comparison of IdPairVec instances"""
        ipv1Copy = ap.IdPairVec(self.ipv1)
        ipv2Copy = ap.IdPairVec(self.ipv2)
        assert ipv1Copy == self.ipv1
        assert ipv2Copy == self.ipv2
        ipv1Copy.swap(ipv2Copy)
        assert ipv1Copy == self.ipv2
        assert ipv2Copy == self.ipv1
        ipv1Copy.swap(ipv2Copy)
        if ipv1Copy.size() == 0:
            ipv1Copy.push_back(ap.IdPair())
        else:
            ipv1Copy.pop_back()
        ip = ap.IdPair()
        ipv2Copy.push_back((123456789876543210, 876543210123456789))
        assert ipv1Copy != self.ipv1
        assert ipv2Copy != self.ipv2

    def testInsertErase(self):
        """Tests inserting and erasing IdPairVec elements"""
        ipv1Copy = ap.IdPairVec(self.ipv1)
        ipv1Copy.insert(ipv1Copy.begin() + 7, ap.IdPair())
        ipv1Copy.insert(ipv1Copy.begin() + 8, 3, ap.IdPair())
        ipv1Copy.erase(ipv1Copy.begin() + 7)
        ipv1Copy.erase(ipv1Copy.begin() + 7, ipv1Copy.begin() + 10)
        assert ipv1Copy == self.ipv1

    def testSlice(self):
        slice = self.ipv1[3:9]
        j = 3
        for i in slice:
            print i
            assert i[0] == j
            j += 1

    def testPersistence(self):
        if persistence.DbAuth.available():
            pol  = policy.PolicyPtr()
            pers = persistence.Persistence.getPersistence(pol)
            loc  = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test")
            dp   = dafBase.DataProperty.createPropertyNode("root")
            dp.addProperty(dafBase.DataProperty("visitId", int(time.clock())*16384 + random.randint(0,16383)))
            dp.addProperty(dafBase.DataProperty("itemName", "IdPair"))
            stl = persistence.StorageList()
            stl.push_back(pers.getPersistStorage("DbStorage", loc))
            pers.persist(self.ipv1, stl, dp)
            stl = persistence.StorageList()
            stl.push_back(pers.getRetrieveStorage("DbStorage", loc))
            res = ap.IdPairVecSharedPtr(pers.retrieve("IdPairVector", stl, dp))
            db = persistence.DbStorage()
            db.setPersistLocation(loc)
            db.startTransaction()
            db.dropTable(ap.getTableName(pol, dp))
            db.endTransaction()
            assert(res == self.ipv1)


# ----------------------------------------------------------------
class IdVecTestCase(unittest.TestCase):
    """A test case for IdVec"""

    def setUp(self):
        n = 16 + random.randint(0,16)
        self.idv1 = ap.IdVec(n)
        self.idv2 = ap.IdVec()

        for i in xrange(n):
            self.idv1[i] = i
            j = (i*17) % 13
            self.idv2.push_back(j)

    def tearDown(self):
        del self.idv1
        del self.idv2

    def testIterable(self):
        """Tests that IdVec instances can be iterated over"""
        j = 0
        for i in self.idv1:
            assert i == j
            j += 1

    def testCopyAndCompare(self):
        """Tests copying, assignment and comparison of IdVec instances"""
        idv1Copy = ap.IdVec(self.idv1)
        idv2Copy = ap.IdVec(self.idv2)
        assert idv1Copy == self.idv1
        assert idv2Copy == self.idv2
        idv1Copy.swap(idv2Copy)
        assert idv1Copy == self.idv2
        assert idv2Copy == self.idv1
        idv1Copy.swap(idv2Copy)
        if idv1Copy.size() == 0:
            idv1Copy.push_back(1409756109)
        else:
            idv1Copy.pop_back()
        j = 876543210123456789
        idv2Copy.push_back(j)
        assert idv1Copy != self.idv1
        assert idv2Copy != self.idv2

    def testInsertErase(self):
        """Tests inserting and erasing IdVec elements"""
        idv1Copy = ap.IdVec(self.idv1)
        idv1Copy.insert(idv1Copy.begin() + 13, 95670927)
        idv1Copy.insert(idv1Copy.begin() + 14, 1, 17465108724650197)
        idv1Copy.erase(idv1Copy.begin() + 13)
        idv1Copy.erase(idv1Copy.begin() + 13, idv1Copy.begin() + 14)
        assert idv1Copy == self.idv1

    def testSlice(self):
        slice = self.idv1[1:11]
        j = 1
        for i in slice:
            print i
            assert i == j
            j += 1

    def testPersistence(self):
        if persistence.DbAuth.available():
            pol  = policy.PolicyPtr()
            pers = persistence.Persistence.getPersistence(pol)
            loc  = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test")
            dp   = dafBase.DataProperty.createPropertyNode("root")
            dp.addProperty(dafBase.DataProperty("visitId", int(time.clock())*16384 + random.randint(0,16383)))
            dp.addProperty(dafBase.DataProperty("itemName", "Id"))
            stl = persistence.StorageList()
            stl.push_back(pers.getPersistStorage("DbStorage", loc))
            pers.persist(self.idv1, stl, dp)
            stl = persistence.StorageList()
            stl.push_back(pers.getRetrieveStorage("DbStorage", loc))
            res = ap.IdVecSharedPtr(pers.retrieve("IdVector", stl, dp))
            db = persistence.DbStorage()
            db.setPersistLocation(loc)
            db.startTransaction()
            db.dropTable(ap.getTableName(pol, dp))
            db.endTransaction()
            assert(res == self.idv1)


# ----------------------------------------------------------------
def suite():
    """Returns a suite containing all the test cases in this module."""

    tests.init()

    suites = []
    suites += unittest.makeSuite(MatchPairVecTestCase)
    suites += unittest.makeSuite(IdPairVecTestCase)
    suites += unittest.makeSuite(IdVecTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    tests.run(suite())

