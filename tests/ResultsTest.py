#! /usr/bin/env python

"""
Tests for C++ association pipeline result vector Python wrappers (including persistence)

Run with:
   python ResultTest.py
or
   python
   >>> import unittest; T=load("ResultTest"); unittest.TextTestRunner(verbosity=1).run(T.suite())
"""

import itertools
import pdb
import random
import time
import unittest

import lsst.daf.base as base
import lsst.pex.policy as policy
import lsst.daf.persistence as persistence
import lsst.afw.detection as detection
import lsst.utils.tests as tests
import lsst.ap as ap

def _seqEqual(seq0, seq1):
    if len(seq0) != len(seq1):
        return False
    for e0, e1 in itertools.izip(seq0, seq1):
        if e0 != e1:
            return False
    return True

def _insertErase(vec, vecType, elemType):
    front = vec[:8]
    back = vec[8:]
    copy = vecType()
    for e in front:
        copy.append(e)
    e = elemType()
    for i in xrange(4):
        copy.append(e)
    for e in back:
        copy.append(e)
    del copy[8]
    del copy[8:11]
    assert _seqEqual(copy, vec)

def _copy(vec, vecType, elemType):
    copy0 = vecType(vec)
    copy1 = vecType(copy0)
    copy1.push_back(elemType())
    assert _seqEqual(copy0, vec)
    assert not _seqEqual(copy1, vec)
    copy0.swap(copy1)
    assert _seqEqual(copy1, vec)
    assert not _seqEqual(copy0, vec)


# ----------------------------------------------------------------
class MatchPairVecTestCase(unittest.TestCase):
    """A test case for MatchPair and MatchPairVec"""

    def setUp(self):
        n = 16 + random.randint(0,16)
        self.vec = ap.MatchPairVec()
        for i in xrange(n):
            mp = ap.MatchPair()
            mp.setFirst(i)
            mp.setSecond(i*20)
            mp.setDistance(random.random())
            self.vec.push_back(mp)

    def tearDown(self):
        del self.vec

    def testIterable(self):
        """Tests that MatchPairVec instances can be iterated over"""
        j = 0
        for i in self.vec:
            assert i.getFirst() == j
            j += 1

    def testCopy(self):
        """Tests copying and assignment of MatchPairVec instances"""
        _copy(self.vec, ap.MatchPairVec, ap.MatchPair)

    def testInsertErase(self):
        """Tests inserting and erasing MatchPairVec elements"""
        _insertErase(self.vec, ap.MatchPairVec, ap.MatchPair)

    def testSlice(self):
        slice = self.vec[1:5]
        j = 1
        for i in slice:
            print i
            assert i.getFirst() == j
            j += 1

    def testPersistence(self):
        if persistence.DbAuth.available("lsst10.ncsa.uiuc.edu", "3306"):
            pol  = policy.Policy()
            root = "Formatter.PersistableMatchPairVector"
            pol.set(root + ".MatchPair.templateTableName", "_tmpl_MatchPair")
            pol.set(root + ".MatchPair.tableNamePattern", "_tmp_v%(visitId)_MP")
            pers = persistence.Persistence.getPersistence(pol)
            loc  = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test")
            ps   = base.PropertySet()
            ps.setInt("visitId", int(time.clock())*16384 + random.randint(0,16383))
            ps.setString("itemName", "MatchPair")
            stl = persistence.StorageList()
            stl.push_back(pers.getPersistStorage("DbStorage", loc))
            pers.persist(ap.PersistableMatchPairVector(self.vec), stl, ps)
            stl = persistence.StorageList()
            stl.push_back(pers.getRetrieveStorage("DbStorage", loc))
            res = ap.PersistableMatchPairVector.swigConvert(pers.unsafeRetrieve("PersistableMatchPairVector", stl, ps))
            db = persistence.DbStorage()
            db.setPersistLocation(loc)
            db.startTransaction()
            db.dropTable(detection.getTableName(pol.getPolicy(root), ps))
            db.endTransaction()
            assert _seqEqual(res.getMatchPairs(), self.vec)


# ----------------------------------------------------------------
class IdPairVecTestCase(unittest.TestCase):
    """A test case for IdPair and IdPairVec"""

    def setUp(self):
        n = 16 + random.randint(0,16)
        self.vec = ap.IdPairVec()
        for i in xrange(n):
            self.vec.push_back((i, -i))

    def tearDown(self):
        del self.vec

    def testIterable(self):
        """Tests that IdPairVec instances can be iterated over"""
        j = 0
        for i in self.vec:
            assert i[0] == j
            j += 1

    def testCopy(self):
        """Tests copying and assignment of IdPairVec instances"""
        _copy(self.vec, ap.IdPairVec, ap.IdPair)

    def testInsertErase(self):
        """Tests inserting and erasing IdPairVec elements"""
        _insertErase(self.vec, ap.IdPairVec, ap.IdPair)

    def testSlice(self):
        slice = self.vec[3:9]
        j = 3
        for i in slice:
            print i
            assert i[0] == j
            j += 1

    def testPersistence(self):
        if persistence.DbAuth.available("lsst10.ncsa.uiuc.edu", "3306"):
            pol  = policy.Policy()
            root = "Formatter.PersistableIdPairVector"
            pol.set(root + ".IdPair.templateTableName", "_tmpl_IdPair")
            pol.set(root + ".IdPair.tableNamePattern", "_tmp_v%(visitId)_IP")
            pers = persistence.Persistence.getPersistence(pol)
            loc  = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test")
            ps   = base.PropertySet()
            ps.setInt("visitId", int(time.clock())*16384 + random.randint(0,16383))
            ps.setString("itemName", "IdPair")
            stl = persistence.StorageList()
            stl.push_back(pers.getPersistStorage("DbStorage", loc))
            pers.persist(ap.PersistableIdPairVector(self.vec), stl, ps)
            stl = persistence.StorageList()
            stl.push_back(pers.getRetrieveStorage("DbStorage", loc))
            res = ap.PersistableIdPairVector.swigConvert(pers.unsafeRetrieve("PersistableIdPairVector", stl, ps))
            db = persistence.DbStorage()
            db.setPersistLocation(loc)
            db.startTransaction()
            db.dropTable(detection.getTableName(pol.getPolicy(root), ps))
            db.endTransaction()
            assert _seqEqual(res.getIdPairs(), self.vec)


# ----------------------------------------------------------------
#class IdVecTestCase(unittest.TestCase):
#    """A test case for IdVec"""
#
#    def setUp(self):
#        n = 16 + random.randint(0,16)
#        self.vec = ap.IdVec()
#        for i in xrange(n):
#            self.vec.push_back(i)
#
#    def tearDown(self):
#        del self.vec
#
#    def testIterable(self):
#        """Tests that IdVec instances can be iterated over"""
#        j = 0
#        for i in self.vec:
#            assert i == j
#            j += 1
#
#    def testCopy(self):
#        """Tests copying and assignment of IdVec instances"""
#        _copy(self.vec, ap.IdVec, long)
#
#    def testInsertErase(self):
#        """Tests inserting and erasing IdVec elements"""
#        _insertErase(self.vec, ap.IdVec, long)
#
#    def testSlice(self):
#        slice = self.vec[1:11]
#        j = 1
#        for i in slice:
#            print i
#            assert i == j
#            j += 1
#
#    def testPersistence(self):
#        if persistence.DbAuth.available("lsst10.ncsa.uiuc.edu", "3306"):
#            pol  = policy.Policy()
#            root = "Formatter.PersistableIdVector"
#            pol.set(root + ".Id.templateTableName", "_tmpl_Id")
#            pol.set(root + ".Id.tableNamePattern", "_tmp_v%(visitId)_I")
#            pers = persistence.Persistence.getPersistence(pol)
#            loc  = persistence.LogicalLocation("mysql://lsst10.ncsa.uiuc.edu:3306/test")
#            ps   = base.PropertySet()
#            ps.setInt("visitId", int(time.clock())*16384 + random.randint(0,16383))
#            ps.setString("itemName", "Id")
#            stl = persistence.StorageList()
#            stl.push_back(pers.getPersistStorage("DbStorage", loc))
#            pers.persist(ap.PersistableIdVector(self.vec), stl, ps)
#            stl = persistence.StorageList()
#            stl.push_back(pers.getRetrieveStorage("DbStorage", loc))
#            res = ap.PersistableIdVector.swigConvert(pers.unsafeRetrieve("PersistableIdVector", stl, ps))
#            db = persistence.DbStorage()
#            db.setPersistLocation(loc)
#            db.startTransaction()
#            db.dropTable(detection.getTableName(pol.getPolicy(root), ps))
#            db.endTransaction()
#            assert _seqEqual(res.getIds(), self.vec)


# ----------------------------------------------------------------
def suite():
    """Returns a suite containing all the test cases in this module."""

    tests.init()

    suites = []
    suites += unittest.makeSuite(MatchPairVecTestCase)
    suites += unittest.makeSuite(IdPairVecTestCase)
    # suites += unittest.makeSuite(IdVecTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    tests.run(suite())

