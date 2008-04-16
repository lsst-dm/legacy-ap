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

import datetime
import os
import os.path
import re
import subprocess

import lsst.pex.exceptions
import lsst.daf.persistence
import lsst.pex.logging
import lsst.pex.harness.Stage
import lsst.afw.image
import lsst.ap as ap


# ----------------------------------------------------------------

class LoadStage(lsst.pex.harness.Stage.Stage):
    """
    Stage that loads basic object data (id, position, variability probabilities)
    for a visit FOV.
    """

    def _massagePolicy(self):
        runDict = { 'runId': str(self.getRun()) }
        if self._policy.exists('filterTableLocation'):
            ftLoc = self._policy.get('filterTableLocation') % runDict
            self._policy.set('filterTableLocation', ftLoc)

    def __init__(self, stageId, policy):
        lsst.pex.harness.Stage.Stage.__init__(self, stageId, policy)
        self.firstVisit = True

    def makeVpContext(self):
        """
        Takes a clipboard from the input queue, creates a visit processing context
        (a holder for all inter-stage association pipeline state), adds it to the
        clipboard, and finally places the clipboard onto the output queue. Expected
        on the input clipboard is an event named 'triggerAssociationEvent' containing
        the following information:
        - the position of the FOV center, given by keys 'FOVRA' and 'FOVDec' (both
          with double precision values in units of degrees).
        - a visit identifier given by the 'visitId' key (with an int64_t value).
        - the name of the filter (a string) for the visit given by the 'filterName' key.
        - [optionally] a match radius given by the 'matchRadius' key (with a double
          precision value in units of arc-seconds).
        """
        clipboard      = self.inputQueue.getNextDataset()
        event          = clipboard.get('triggerAssociationEvent')
        self.vpContext = ap.VisitProcessingContext(
            event, self.getRun(), self.getRank(), self.getUniverseSize() - 1
        )
        clipboard.put('vpContext', self.vpContext)
        self.outputQueue.addDataset(clipboard)

    def preprocess(self):
        """
        Registers the incoming visit with the shared memory chunk manager.
        """
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        lsst.pex.logging.Trace('associate.LoadStage', 3, 'Python lsst.ap.pipeline.LoadStage preprocess(): stage %d' % self.getStageId())

        if self.firstVisit:
            self._massagePolicy()
            ap.initialize(self._policy, str(self.getRun()))
            self.firstVisit = False
        self.makeVpContext()
        ap.registerVisit(self.vpContext)

    def process(self):
        """
        Loads object chunk and chunk delta files assigned to the slice.
        """
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        lsst.pex.logging.Trace('associate.LoadStage', 3, 'Python lsst.ap.pipeline.LoadStage process(): stage %d' % self.getStageId())
        lsst.pex.logging.Trace('associate.LoadStage', 3, 'Python lsst.ap.pipeline.LoadStage process(): worker %d' % self.getRank())

        if self.firstVisit:
            self._massagePolicy()
            ap.initialize(self._policy, str(self.getRun()))
            self.firstVisit = False
        self.makeVpContext()
        ap.loadSliceObjects(self.vpContext)

    def postprocess(self):
        """
        Checks to make sure all worker slices successfully loaded their share of the
        objects for the visit FOV and builds an object index if so.
        """
        lsst.pex.logging.Trace('associate.LoadStage', 3, 'Python lsst.ap.pipeline.LoadStage postprocess(): stage %d' % self.getStageId())
        ap.buildObjectIndex(self.vpContext)


# --------------------------------------------------------------------------------

class MatchDiaSourcesStage(lsst.pex.harness.Stage.Stage):
    """
    Matches difference sources from the detection pipeline against objects
    within the visits FOV. Difference sources are expected to be found on the clipboard
    under the key 'DiaSources' and match results, consisting of (difference source id,
    object id, match distance) tuples are placed onto the clipboard under the key
    'DiaSourceToObjectMatches'.
    """

    def __init__(self, stageId=-1, policy=None):
        lsst.pex.harness.Stage.Stage.__init__(self, stageId, policy)

    def preprocess(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        lsst.pex.logging.Trace('associate.MatchDiaSourcesStage', 3, 'Python lsst.ap.pipeline.MatchDiaSourcesStage preprocess(): stage %d' % self.getStageId())

        clipboard = self.inputQueue.getNextDataset()
        vpContext = clipboard.get('vpContext')
        sources   = clipboard.get('DiaSources')
        matches   = ap.MatchPairVecPtr()

        vpContext.setDiaSources(sources)
        ap.matchDiaSources(matches, vpContext)

        clipboard.put('DiaSourceToObjectMatches', matches)
        self.outputQueue.addDataset(clipboard)

    def process(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0

        self.outputQueue.addDataset(self.inputQueue.getNextDataset())

    def postprocess(self):
        pass


# --------------------------------------------------------------------------------

class MatchMopsPredsStage(lsst.pex.harness.Stage.Stage):
    """
    Matches moving object predictions for a visit against difference sources. The previous
    stage is assumed to have read the difference sources and produced an index for them (stored
    in 'vpContext' on the clipboard). Moving object predictions are expected to be found on the
    clipboard under the key 'MopsPreds' and match results, consisting of (moving object id,
    difference source id, match distance) tuples are placed onto the clipboard under the key
    'MopsPredToDiaSourceMatches'. Finally, difference sources which didn't match anything are
    used to create new objects - identifiers for these are placed onto the clipboard under the
    key 'NewObjects'.
    """
    def __init__(self, stageId=-1, policy=None):
        lsst.pex.harness.Stage.Stage.__init__(self, stageId, policy)

    def preprocess(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        lsst.pex.logging.Trace('associate.MatchMopsPredsStage', 3, 'Python lsst.ap.pipeline.MatchMopsPredsStage preprocess(): stage %d' % self.getStageId())

        clipboard = self.inputQueue.getNextDataset()
        vpContext = clipboard.get('vpContext')

        try:
            event1 = clipboard.get('triggerAssociationEvent')
            event2 = clipboard.get('triggerMatchMopsPredsEvent')
            if event1.findUnique('visitId').getValueInt() != event2.findUnique('visitId').getValueInt():
                raise lsst.pex.exceptions.LsstRuntime(
                    'triggerAssociationEvent.visitId != triggerMatchMopsPredsEvent.visitId'
                )
            preds     = clipboard.get('MopsPreds')
            matches   = ap.MatchPairVecPtr()
            idPairs   = ap.IdPairVecPtr()
        except:
            ap.endVisit(vpContext, True)
            raise

        ap.matchMops(matches, idPairs, vpContext, preds)

        clipboard.put('MopsPredToDiaSourceMatches', matches)
        clipboard.put('NewObjectIdPairs', idPairs)
        self.outputQueue.addDataset(clipboard)

    def process(self):
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0

        self.outputQueue.addDataset(self.inputQueue.getNextDataset())

    def postprocess(self):
        pass


# ----------------------------------------------------------------

class StoreStage(lsst.pex.harness.Stage.Stage):
    """
    Stores new objects (created from unmatched difference sources) obtained during
    the visit into chunk delta files.
    """

    def _createSqlScripts(self, vpContext, visitTime):
        self.templateDict['visitId']    = vpContext.getVisitId()
        self.templateDict['visitTime']  = visitTime
        self.templateDict['filterName'] = self.filterChars[vpContext.getFilterId()]
        for i in self.templateNames:
            outText = self.templates[i] % self.templateDict
            name, ext = os.path.splitext(i)
            assert name.endswith('Template')
            outName = ''.join([name[:-len('Template')], '_visit%d' % vpContext.getVisitId(), ext])
            outPath = os.path.join(self.scriptDir, outName)
            self.scriptPaths[i] = outPath
            with file(outPath, 'w') as outFile:
                outFile.write(outText)

    def _runSqlScript(self, i):
        script = self.scriptPaths[i]
        args = ['mysql',
                '--user=%s' % self.username,
                '--password=%s' % self.password,
                '--host=%s' % self.hostname,
                '--port=%d' % self.port,
                '--database=%s' % self.database,
                '-vvv',
                '-e', 'source %s' % script]
        with file(script + '.log', 'w') as logFile:
            exitcode = subprocess.Popen(args, stdout=logFile).wait()
        return exitcode

    def _initPropsFromRun(self, runId):
        runDict  = { 'runId': str(runId) }
        location = self.location % runDict
        self.scriptDir = self.scriptDirectory % runDict
        if not os.path.exists(self.scriptDir):
            os.makedirs(self.scriptDir)
        # Parse database location string
        dbloc = re.compile('(\w+)://(\S+):(\d+)/(\S+)').match(location)
        if dbloc is None:
            raise lsst.pex.exceptions.LsstRuntime('invalid location string')
        if dbloc.group(1) != 'mysql':
            raise lsst.pex.exceptions.LsstRuntime('database type %s not supported' % dbloc.group(1))
        self.hostname = dbloc.group(2)
        self.port     = int(dbloc.group(3))
        self.database = dbloc.group(4)

    def __init__(self, stageId=-1, policy=None):
        lsst.pex.harness.Stage.Stage.__init__(self, stageId, policy)
        if not lsst.daf.persistence.DbAuth.available():
            raise lsst.pex.exceptions.LsstRuntime('missing credentials for database authorization')
        self.firstVisit    = True
        self.username      = lsst.daf.persistence.DbAuth.username()
        self.password      = lsst.daf.persistence.DbAuth.password()
        self.filterChars   = ('u','g','r','i','z','y')
        self.templateNames = ['StoreOutputsTemplate.sql',
                              'AppendTablesTemplate.sql',
                              'DropTablesTemplate.sql']
        self.templates     = {}
        self.scriptPaths   = {}
        self.templateDict  = {}
        # read in templates
        associateDir = os.environ['ASSOCIATE_DIR']
        for i in self.templateNames:
            inPath = os.path.join(associateDir, "sql", i)
            with file(inPath, 'r') as inFile:
                self.templates[i] = inFile.read()
        # set default policy parameters
        self.location        = 'mysql://lsst10.ncsa.uiuc.edu:3306/test'
        self.storeOutputs    = True
        self.appendTables    = True
        self.dropTables      = True
        self.scriptDirectory = '/tmp/%(runId)s'
        self.templateDict['diaSourceTable']    = 'DIASource'
        self.templateDict['varObjectTable']    = 'VarObject'
        self.templateDict['nonVarObjectTable'] = 'NonVarObject'
        # extract policy parameters from policy (if available)
        if policy != None:
            self.location        = policy.getString('location', self.location)
            self.storeOutputs    = policy.getBool('storeOutputs', self.storeOutputs)
            self.appendTables    = policy.getBool('appendTables', self.appendTables)
            self.dropTables      = policy.getBool('dropTables', self.dropTables)
            self.scriptDirectory = policy.getString('scriptDirectory', self.scriptDirectory)
            self.templateDict['diaSourceTable'] = policy.getString('diaSourceTable', self.templateDict['diaSourceTable'])
            self.templateDict['varObjectTable'] = policy.getString('varObjectTable', self.templateDict['varObjectTable'])
            self.templateDict['nonVarObjectTable'] = policy.getString('nonVarObjectTable', self.templateDict['nonVarObjectTable'])
        self._initPropsFromRun(self.getRun())

    def preprocess(self):
        if self.firstVisit:
            self._initPropsFromRun(self.getRun())
            self.firstVisit = False

    def process(self):
        """
        Stores chunk deltas assigned to the worker slice.
        """
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        lsst.pex.logging.Trace('associate.StoreStage', 3, 'Python lsst.ap.pipeline.StoreStage process(): stage %d' % self.getStageId())
        lsst.pex.logging.Trace('associate.StoreStage', 3, 'Python lsst.ap.pipeline.StoreStage process(): worker %d' % self.getRank())

        clipboard = self.inputQueue.getNextDataset()
        self.outputQueue.addDataset(clipboard)
        ap.storeSliceObjects(clipboard.get('vpContext'))

    def postprocess(self):
        """
        Checks to make sure all worker slices successfully stored their share of the
        new objects for the visit and removes the visit from the list of in-flight
        visits being tracked by the shared memory chunk manager. Optionally:
          - stores pipeline outputs (new objects, object updates, new difference sources)
          - appends per-visit tables to global accumulator tables
          - drops per-visit tables created by LSST pipelines
        """
        assert self.inputQueue.size()  == 1
        assert self.outputQueue.size() == 0
        lsst.pex.logging.Trace('associate.StoreStage', 3, 'Python lsst.ap.pipeline.StoreStage postprocess(): stage %d' % self.getStageId())

        clipboard = self.inputQueue.getNextDataset()
        self.outputQueue.addDataset(clipboard)
        vpContext = clipboard.get('vpContext')
        event     = clipboard.get('triggerAssociationEvent')
        # get MJD of visit and convert to UTC string in ISO 8601 format for use in database queries
        dt        = lsst.daf.persistence.DateTime(event.findUnique('visitTime').getValueDouble())
        utcString = datetime.datetime.utcfromtimestamp(dt.nsecs()/1000000000).isoformat(' ')
        self._createSqlScripts(vpContext, utcString)
        if self.storeOutputs:
            if self._runSqlScript('StoreOutputsTemplate.sql') != 0:
                ap.endVisit(vpContext, True)
                raise lsst.pex.exceptions.LsstRuntime(
                    'Association pipeline failed: SQL script %s failed' %
                    self.scriptPaths['StoreOutputsTemplate.sql']
                )

        if not ap.endVisit(vpContext, False):
            raise lsst.pex.exceptions.LsstRuntime('Association pipeline failed: visit not committed')

        if self.appendTables:
            if self._runSqlScript('AppendTablesTemplate.sql') != 0:
                raise lsst.pex.exceptions.LsstRuntime(
                    'Association pipeline failed: SQL script %s failed' %
                    self.scriptPaths['AppendTablesTemplate.sql']
                )

        if self.dropTables:
            if self._runSqlScript('DropTablesTemplate.sql') != 0:
                raise lsst.pex.exceptions.LsstRuntime(
                    'Association pipeline failed to drop tables given in SQL script %s' %
                    self.scriptPaths['DropTablesTemplate.sql']
                )


# ----------------------------------------------------------------

