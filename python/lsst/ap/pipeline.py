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

"""
Classes corresponding to the 4 Stages that make up the LSST association pipeline:

 - load pre-existing object data for the visit FOV from chunk files
 - match difference sources for the visit (from the detection pipeline) against objects,
   outputting match results to the database
 - match moving object predictions for the visit (from MOPS) against difference sources,
   outputting match results and new objects (created from difference sources with no matches)
   to the database
 - store new objects into chunk delta files
"""
from __future__ import with_statement

from datetime import datetime
import os, os.path
import pdb
import re
import time

import lsst.daf.base as base
import lsst.daf.persistence as persistence
import lsst.pex.logging as logging
import lsst.pex.harness as harness
import lsst.afw.image as image
import lsst.ap as ap


# ----------------------------------------------------------------

class LoadStage(harness.Stage.Stage):
    """
    Stage that loads basic object data (id, position, variability probabilities)
    for a visit FOV.
    """

    def _massagePolicy(self):
        loc = persistence.LogicalLocation(self._policy.get('filterTableLocation'))
        self._policy.set('filterTableLocation', loc.locString())

    def __init__(self, stageId, policy):
        if policy is None:
            raise RuntimeError, "Cannot create a lsst.ap.LoadStage without a policy"
        harness.Stage.Stage.__init__(self, stageId, policy)
        self._firstVisit = True

    def makeVpContext(self):
        """
        Takes a clipboard from the input queue, creates a visit processing context
        (a holder for all inter-stage association pipeline state), adds it to the
        clipboard, and finally places the clipboard onto the output queue. Expected
        on the input clipboard is an event (named 'triggerAssociationEvent')
        containing the following information:
        - the position of the visit center, given by keys 'ra' and 'decl' (both
          with double precision values in units of degrees).
        - a visit identifier given by the 'visitId' key (with an integer value).
        - the name of the filter (a string) for the visit given by the 'filter' key.
        - the time at which the visit will occur, given by key 'dateObs' with a double
          precision value in MJD (TAI)
        - [optionally] a match radius given by the 'matchRadius' key (with a double
          precision value in units of arc-seconds).
        """
        clipboard = self.inputQueue.getNextDataset()
        event = clipboard.get('triggerAssociationEvent')
        self.vpContext = ap.VisitProcessingContext(
            self._policy, event, self.getRun(), self.getRank(), self.getUniverseSize() - 1
        )
        clipboard.put('vpContext', self.vpContext)
        self.outputQueue.addDataset(clipboard)

    def preprocess(self):
        """
        Registers the incoming visit with the shared memory chunk manager.
        """
        assert self.inputQueue.size() == 1
        assert self.outputQueue.size() == 0
        if self._firstVisit:
            self._massagePolicy()
            ap.initialize(str(self.getRun()))
            self._firstVisit = False
        self.makeVpContext()
        ap.registerVisit(self.vpContext)

    def process(self):
        """
        Loads object chunk and chunk delta files assigned to the slice.
        """
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        if self._firstVisit:
            self._massagePolicy()
            ap.initialize(str(self.getRun()))
            self._firstVisit = False
        self.makeVpContext()
        ap.loadSliceObjects(self.vpContext)

    def postprocess(self):
        """
        Checks to make sure all worker slices successfully loaded their share of the
        objects for the visit FOV and builds an object index if so.
        """
        ap.buildObjectIndex(self.vpContext)


# --------------------------------------------------------------------------------

class MatchDiaSourcesStage(harness.Stage.Stage):
    """
    Matches difference sources from the detection pipeline against objects
    within the visits FOV. Difference sources are expected to be found on the clipboard
    under the key 'diaSources' and match results, consisting of (difference source id,
    object id, match distance) tuples are placed onto the clipboard under the key
    'diaSourceToObjectMatches'.
    """

    def __init__(self, stageId=-1, policy=None):
        harness.Stage.Stage.__init__(self, stageId, policy)

    def preprocess(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        clipboard = self.inputQueue.getNextDataset()
        vpContext = clipboard.get('vpContext')
        sources   = clipboard.get('diaSources')
        matches   = ap.MatchPairVec()

        vpContext.setDiaSources(sources)
        ap.matchDiaSources(matches, vpContext)

        clipboard.put('diaSourceToObjectMatches', ap.PersistableMatchPairVector(matches))
        self.outputQueue.addDataset(clipboard)

    def process(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        self.outputQueue.addDataset(self.inputQueue.getNextDataset())

    def postprocess(self):
        pass


# --------------------------------------------------------------------------------

class MatchMopsPredsStage(harness.Stage.Stage):
    """
    Matches moving object predictions for a visit against difference sources. The previous
    stage is assumed to have read the difference sources and produced an index for them (stored
    in 'vpContext' on the clipboard). Moving object predictions are expected to be found on the
    clipboard under the key 'mopsPreds' and match results, consisting of (moving object id,
    difference source id, match distance) tuples are placed onto the clipboard under the key
    'predToDiaSourceMatches'. Finally, difference sources which didn't match anything are
    used to create new objects - identifiers for these are placed onto the clipboard under the
    key 'diaSourceToNewObject'.
    """
    def __init__(self, stageId=-1, policy=None):
        harness.Stage.Stage.__init__(self, stageId, policy)

    def preprocess(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        clipboard = self.inputQueue.getNextDataset()
        vpContext = clipboard.get('vpContext')

        try:
            ev1 = clipboard.get('triggerAssociationEvent')
            ev2 = clipboard.get('triggerMatchMopsPredsEvent')
            if ev1.getInt('visitId') != ev2.getInt('visitId'):
                raise RuntimeError('triggerAssociationEvent.visitId != triggerMatchMopsPredsEvent.visitId')
            preds = clipboard.get('mopsPreds')
            matches = ap.MatchPairVec()
            idPairs = ap.IdPairVec()
        except:
            ap.endVisit(vpContext, True)
            raise

        ap.matchMops(matches, idPairs, vpContext, preds.getPredictions())

        clipboard.put('predToDiaSourceMatches', ap.PersistableMatchPairVector(matches))
        clipboard.put('diaSourceToNewObject', ap.PersistableIdPairVector(idPairs))
        self.outputQueue.addDataset(clipboard)

    def process(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        self.outputQueue.addDataset(self.inputQueue.getNextDataset())

    def postprocess(self):
        pass


# ----------------------------------------------------------------

class StoreStage(harness.Stage.Stage):
    """
    Store new objects (created from unmatched difference sources under certain
    conditions) obtained during the visit into chunk delta files.
    """
    def _runSql(self, sqlStatements, scriptFileName):
        db = persistence.DbStorage()
        db.setPersistLocation(self.database)
        with open(os.path.join(self.scriptDir, scriptFileName), 'w') as scriptFile:
            for stmt in sqlStatements:
                stmt = stmt % self.templateDict
                startTime = time.clock()
                scriptFile.write(stmt)
                scriptFile.write(';\n\n')
                scriptFile.flush()
                try:
                    db.executeSql(stmt)
                except Exception, e:
                    try:
                        for line in str(e).splitlines(True):
                            scriptFile.write('-- ' + line)
                        scriptFile.flush()
                    except:
                        pass
                    raise
                scriptFile.write('-- statement executed in %f seconds\n\n' % (time.clock() - startTime))

    def _copyObject(self):
        """
        Append contents of policy-specified object table to Object table in per-run database
        """
        if (self.objectTable and len(self.objectTable) > 0):
            db = persistence.DbStorage()
            db.setPersistLocation(self.database)
            db.executeSql("INSERT INTO %s SELECT * FROM %s" %\
                (self.templateDict['nonVarObjectTable'], self.objectTable))

    def __init__(self, stageId, policy):
        if policy is None:
            raise RuntimeError, "Cannot create a lsst.ap.StoreStage without a policy"
        harness.Stage.Stage.__init__(self, stageId, policy)
        self.filterChars  = ('u','g','r','i','z','y')
        self.templateDict = {}
        self.additionalData = base.PropertySet()
        self.database = persistence.LogicalLocation(policy.getString('database'))
        self.objectTable = policy.getString('objectTable') if policy.exists('objectTable') else None
        self.storeOutputs = policy.getBool('storeOutputs')
        self.appendTables = policy.getBool('appendTables')
        self.dropTables = policy.getBool('dropTables')
        self.scriptDir = persistence.LogicalLocation(policy.getString('scriptDirectory')).locString()
        if not os.path.exists(self.scriptDir):
            os.makedirs(self.scriptDir)
        self.templateDict['diaSourceTable'] = policy.getString('diaSourceTable')
        self.templateDict['varObjectTable'] = policy.getString('varObjectTable')
        self.templateDict['nonVarObjectTable'] = policy.getString('nonVarObjectTable')

    def initialize(self, outQueue, inQueue):
        """
        The master slice is in charge of database prep-work
        """
        harness.Stage.Stage.initialize(self, outQueue, inQueue)
        if self.getRank() == -1:
            self._copyObject()

    def preprocess(self):
        pass

    def process(self):
        """
        Store chunk deltas assigned to the worker slice.
        """
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        clipboard = self.inputQueue.getNextDataset()
        self.outputQueue.addDataset(clipboard)
        ap.storeSliceObjects(clipboard.get('vpContext'))

    def postprocess(self):
        """
        Check to make sure all worker slices successfully stored their share of the
        new objects for the visit and removes the visit from the list of in-flight
        visits being tracked by the shared memory chunk manager. Optionally:
          - store pipeline outputs (new objects, object updates, new difference sources)
          - append per-visit tables to global accumulator tables
          - drop per-visit tables created by LSST pipelines
        """
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        clipboard = self.inputQueue.getNextDataset()
        self.outputQueue.addDataset(clipboard)

        vpContext = clipboard.get('vpContext')
        event = clipboard.get('triggerAssociationEvent')
        do = base.DateTime(event.getDouble('dateObs'))
        doUtc = datetime.utcfromtimestamp(do.nsecs(base.DateTime.UTC)/1000000000).isoformat(' ')
        dbType = re.match(r'(\w+)://', self.database.locString()).group(1)
        visitId = event.getInt('visitId')

        self.templateDict['runId'] = self.getRun()
        self.templateDict['visitId'] = visitId
        self.templateDict['dateObs'] = doUtc
        self.templateDict['filter'] = self.filterChars[vpContext.getFilterId()]

        if self.storeOutputs:
            try:
                self._runSql(ap.SqlStoreOutputs.sqlStatements[dbType], 'StoreOutputs_%d.sql' % visitId)
            except:
                ap.endVisit(vpContext, True)
                raise
        if not ap.endVisit(vpContext, False):
            raise RuntimeError('Association pipeline failed: visit not committed')
        if self.appendTables:
            self._runSql(ap.SqlAppendTables.sqlStatements[dbType], 'AppendTables_%d.sql' % visitId)
        if self.dropTables:
            self._runSql(ap.SqlDropTables.sqlStatements[dbType], 'DropTables_%d.sql' % visitId)


# ----------------------------------------------------------------

