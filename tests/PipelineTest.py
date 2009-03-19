#! /usr/bin/env python

"""
Tests a run of the association pipeline on a single visit. This test assumes the existence
of appropriate input tables in the 'test_ap' database. It is intended as a convenient way to
test/debug AP stages without having to deal with MPI and events.

Run using `python PipelineTest.py` or as part of the association pipeline unit tests.

Note that one can arrange for the shared memory object used by the chunk manager to survive the
death of the pipeline process(es). The ap_shmem_admin tool in the bin/ subdirectory of the association
pipeline can then be used to query the chunk manager (--help for details). For this kind of debugging,
the line reading:

    SharedObjectChunkManager::destroyInstance(context.getRunId());

inside the

    void registerVisit(VisitProcessingContext &)

function definition in src/Stages.cc must be commented out. This line normally causes the shared memory
object to be unlinked as soon as all slices have obtained a reference to it and before processing of
the first visit is started, meaning that ap_shmem_admin invocations will only display information on
freshly created (therefore empty) shared memory objects. To remove a left-over shared memory object
one can run

    ap_shmem_admin -u
"""
from itertools import izip
import os, os.path
import pdb
import sys
import time
import unittest

from lsst.daf.base import cout, Citizen, DateTime, PropertySet
from lsst.pex.policy import Policy
from lsst.daf.persistence import DbAuth, DbStorage, LogicalLocation
from lsst.pex.harness.Clipboard import Clipboard
from lsst.pex.harness.Queue import Queue
from lsst.pex.harness.IOStage import InputStage, OutputStage

import lsst.pex.logging as log
import lsst.ap.pipeline as ap


class PipelineTestCase(unittest.TestCase):
    """
    Runs a single visit through the association pipeline. Assumes the existence
    of a test_ap database that contains copies of the necessary input data
    """
    def setUp(self):
        # Turn on tracing
        log.Trace.setVerbosity('', 10)
        log.ScreenLog.createDefaultLog(True, log.Log.INFO)

        # Eventually, these should be read from a policy somewhere
        self.dbServer = 'lsst10.ncsa.uiuc.edu'
        self.dbPort = '3306'
        self.dbType = 'mysql'
        if not DbAuth.available(self.dbServer, self.dbPort):
            self.fail("Cannot access database server %s:%s" % (self.dbServer, self.dbPort))
        # Construct test run database name
        self.runId = DbAuth.username(self.dbServer, self.dbPort) +\
                     time.strftime("_test_ap_%y%m%d_%H%M%S", time.gmtime())

        # Tweak these to run on different input data, or with a different number of slices
        self.universeSize = 2
        self.visitId = 708125
        self.filter = 'u'
        self.ra = 333.880166667
        self.dec = -17.7374166667

        self.dbUrlPrefix = ''.join([self.dbType, '://', self.dbServer, ':', self.dbPort, '/'])
        self.dbUrl = self.dbUrlPrefix + self.runId
        self.substitutions = { 'visitId': self.visitId,
                               'filter': self.filter,
                               'runId': self.runId }
        # Create a database specifically for the test (copy relevant
        # tables from the test_ap database)
        mysqlStatements = [
            """CREATE DATABASE %(runId)s""",
            """USE %(runId)s""",
            """CREATE TABLE Object LIKE test_ap.Object""",
            """CREATE TABLE NonVarObject LIKE test_ap.Object""",
            """CREATE TABLE DIASource LIKE test_ap.DIASource""",
            """CREATE TABLE prv_Filter LIKE test_ap.prv_Filter""",
            """INSERT INTO prv_Filter SELECT * FROM test_ap.prv_Filter""",
            """CREATE TABLE _tmp_v%(visitId)d_DIASource
               LIKE test_ap._tmp_v%(visitId)d_DIASource""",
            """INSERT INTO _tmp_v%(visitId)d_DIASource
               SELECT * FROM test_ap._tmp_v%(visitId)d_DIASource""",
            """CREATE TABLE _tmp_v%(visitId)d_Preds
               LIKE test_ap._tmp_v%(visitId)d_Preds""",
            """INSERT INTO _tmp_v%(visitId)d_Preds
               SELECT * FROM test_ap._tmp_v%(visitId)d_Preds""",
            """CREATE TABLE _tmpl_MatchPair LIKE test_ap._tmpl_MatchPair""",
            """CREATE TABLE _tmpl_IdPair LIKE test_ap._tmpl_IdPair""",
            """CREATE TABLE _tmpl_InMemoryObject LIKE test_ap._tmpl_InMemoryObject""",
            """CREATE TABLE _tmpl_InMemoryMatchPair LIKE test_ap._tmpl_InMemoryMatchPair""",
            """CREATE TABLE _tmpl_InMemoryId LIKE test_ap._tmpl_InMemoryId""",
            """CREATE TABLE _ap_DIASourceToObjectMatches LIKE test_ap._ap_DIASourceToObjectMatches""",
            """CREATE TABLE _ap_PredToDIASourceMatches LIKE test_ap._ap_PredToDIASourceMatches""",
            """CREATE TABLE _ap_DIASourceToNewObject LIKE test_ap._ap_DIASourceToNewObject"""
        ]
        db = DbStorage()
        db.setPersistLocation(LogicalLocation(self.dbUrlPrefix + 'test_ap'))
        try:
            for stmt in mysqlStatements:
                db.executeSql(stmt % self.substitutions)
            
            # Specify list of stages ...
            self.stages = [ ap.LoadStage,
                            InputStage,
                            ap.MatchDiaSourcesStage,
                            OutputStage,
                            InputStage,
                            ap.MatchMopsPredsStage,
                            OutputStage,
                            ap.StoreStage ]

            # and read in stage policy for each stage
            policyDir = os.path.join(os.environ['AP_DIR'], 'pipeline', 'examples', 'policy')
            self.policies = [ Policy(os.path.join(policyDir,'LoadStage.paf')),
                              Policy(os.path.join(policyDir,'MatchDiaSourcesStageInput.paf')),
                              None,
                              Policy(os.path.join(policyDir,'MatchDiaSourcesStageOutput.paf')),
                              Policy(os.path.join(policyDir,'MatchMopsPredsStageInput.paf')),
                              None,
                              Policy(os.path.join(policyDir,'MatchMopsPredsStageOutput.paf')),
                              Policy(os.path.join(policyDir,'StoreStage.paf')) ]

            # construct PropertySet for string interpolation
            psSubs = PropertySet()
            psSubs.setInt('visitId', self.visitId)
            psSubs.setString('runId', self.runId)
            psSubs.setString('filter', self.filter)
            psSubs.setString('work', '.')
            psSubs.setString('input', '/tmp')
            psSubs.setString('output', '/tmp')
            psSubs.setString('update', '/tmp')
            psSubs.setString('dbUrl', self.dbUrl)
            LogicalLocation.setLocationMap(psSubs)
        except:
            # cleanup database in case of error
            db.executeSql("DROP DATABASE %(runId)s" % self.substitutions)
            raise

    def testOneVisit(self):
        # Create a list of clipboards, stage lists, and queue lists for each slice
        clipboards = [Clipboard() for i in xrange(self.universeSize)]
        stageLists = [[] for i in xrange(self.universeSize)]
        queueLists = []
        for i in xrange(self.universeSize):
            queueList = [Queue() for j in xrange(len(self.stages) + 1)]
            queueList[0].addDataset(clipboards[i])
            queueLists.append(queueList)

        # Create and initialize stages for each slice
        for stageClass, policy, i in izip(self.stages, self.policies, xrange(len(self.stages))):
            for stageList, queueList, rank in izip(stageLists, queueLists, xrange(self.universeSize)):
                stage = stageClass(i, policy)
                stage.setRun(self.runId)
                stage.setUniverseSize(self.universeSize)
                stage.setRank(rank - 1)
                stage.initialize(queueList[i+1], queueList[i])
                stageList.append(stage)

        # Create the association pipeline trigger event
        dateObs = DateTime.now().mjd(DateTime.TAI)
        triggerAssociationEvent = PropertySet()
        triggerAssociationEvent.setInt('visitId', self.visitId)
        triggerAssociationEvent.setDouble('dateObs', dateObs)
        triggerAssociationEvent.setString('filter', self.filter)
        triggerAssociationEvent.setDouble('ra', self.ra)
        triggerAssociationEvent.setDouble('decl', self.dec)

        # Create the event triggering the match against moving object predictions
        triggerMatchMopsPredsEvent = PropertySet()
        triggerMatchMopsPredsEvent.setInt('visitId', self.visitId)

        # Add the events to clipboard of each stage
        for clip in clipboards:
            clip.put('triggerAssociationEvent', triggerAssociationEvent)
            clip.put('triggerMatchMopsPredsEvent', triggerMatchMopsPredsEvent)

        assert self.universeSize > 1
        masterStageList = stageLists.pop(0)

        # Run the pipeline (worker slices are run one after the other)
        for masterStage, workerStages in izip(masterStageList, izip(*stageLists)):
            masterStage.preprocess()
            map(lambda x: x.process(), workerStages)
            masterStage.postprocess()

        # Close log to avoid bogus memory-leak reports
        log.Log.closeDefaultLog()

    def tearDown(self):
        """Clean up after test case runs"""
        db = DbStorage()
        db.setPersistLocation(LogicalLocation(self.dbUrlPrefix + 'test_ap'))
        #db.executeSql("DROP DATABASE %(runId)s" % self.substitutions)
        del self.policies
        del self.stages


class MemoryTestCase(unittest.TestCase):
    """Check for memory leaks of citizens"""
    def testLeak(self):
        nleak = Citizen.census(0, 0)
        if nleak != 0:
             Citizen.census(cout, 0)
             self.fail("Leaked %d blocks" % nleak)


def suite():
    """Returns a suite containing all the test cases in this module."""
    suites = [ unittest.makeSuite(PipelineTestCase),
               unittest.makeSuite(MemoryTestCase) ]
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    status = 0 if unittest.TextTestRunner().run(suite()).wasSuccessful() else 1
    return sys.exit(status) if exit else status

if __name__ == "__main__":
    run(True)

