#! /usr/bin/env python

"""
Tests a run of the association pipeline on a single visit. This test enforces a visit id of 1
and assumes the existence of corresponding input tables. The test is therefore not suitable
as a unit test that automatically runs through scons (multiple users running the test concurrently
will interfere with eachother). It is intended as a convenient way to test/debug AP stages without
having to deal with MPI and events.

Run with:
   python AssociationPipelineTest.py

If something is particularly badly broken, some manual cleanup might be required. In particular, some
database tables may need to be dropped and some truncated. Run:
    python TestCleanup.py
to do this.

Note also that one can arrange for the shared memory object used by the chunk manager to survive the
death of the pipeline process(es). The ap_shmem_admin tool in the bin/ subdirectory of the association
pipeline can then be used to query the chunk manager (--help for details). For this kind of debugging,
the line reading:

    SharedSimpleObjectChunkManager::destroyInstance(context.getRunId());

inside the

    void registerVisit(VisitProcessingContext &)

function definition in src/Stages.cc must be commented out. This line normally causes the shared memory
object to be unlinked as soon as all slices have obtained a reference to it and before processing of
the first visit is started, meaning that ap_shmem_admin invocations will only display information on
freshly created (therefore empty) shared memory objects. To remove a left-over shared memory object
one can run

    ap_shmem_admin -u
"""
import os
import os.path
import pdb
import sys
import time
import lsst.dps.Clipboard
import lsst.dps.Queue
import lsst.dps.IOStage
import lsst.mwi.data
import lsst.mwi.exceptions
import lsst.mwi.persistence
import lsst.mwi.policy
import lsst.mwi.utils
import lsst.ap.interface
import lsst.ap.pipeline
import TestCleanup


def runOneVisit():

    runId = 'test'

    # Turn on tracing
    lsst.mwi.utils.Trace.setVerbosity('', 10)
    lsst.mwi.logging.ScreenLog.createDefaultLog(True, lsst.mwi.logging.Log.INFO)

    # Read in stage policies
    policyDir          = os.path.join(os.environ['ASSOCIATE_DIR'], 'pipeline', 'examples', 'policy')
    loadPolicy         = lsst.mwi.policy.Policy.createPolicy(os.path.join(policyDir,'LoadStage.paf'))
    match1InputPolicy  = lsst.mwi.policy.Policy.createPolicy(os.path.join(policyDir,'MatchDiaSourcesStageInput.paf'))
    match1OutputPolicy = lsst.mwi.policy.Policy.createPolicy(os.path.join(policyDir,'MatchDiaSourcesStageOutput.paf'))
    match2InputPolicy  = lsst.mwi.policy.Policy.createPolicy(os.path.join(policyDir,'MatchMopsPredsStageInput.paf'))
    match2OutputPolicy = lsst.mwi.policy.Policy.createPolicy(os.path.join(policyDir,'MatchMopsPredsStageOutput.paf'))
    storePolicy        = lsst.mwi.policy.Policy.createPolicy(os.path.join(policyDir,'StoreStage.paf'))

    # Create stage queues (worker slices and the pipeline master get seperate queues)
    mq0 = lsst.dps.Queue.Queue()
    mq1 = lsst.dps.Queue.Queue()
    mq2 = lsst.dps.Queue.Queue()
    mq3 = lsst.dps.Queue.Queue()
    mq4 = lsst.dps.Queue.Queue()
    mq5 = lsst.dps.Queue.Queue()
    mq6 = lsst.dps.Queue.Queue()
    mq7 = lsst.dps.Queue.Queue()
    mq8 = lsst.dps.Queue.Queue()

    wq0 = lsst.dps.Queue.Queue()
    wq1 = lsst.dps.Queue.Queue()
    wq2 = lsst.dps.Queue.Queue()
    wq3 = lsst.dps.Queue.Queue()
    wq4 = lsst.dps.Queue.Queue()
    wq5 = lsst.dps.Queue.Queue()
    wq6 = lsst.dps.Queue.Queue()
    wq7 = lsst.dps.Queue.Queue()
    wq8 = lsst.dps.Queue.Queue()

    # Create and initialize association pipeline stages - simulate 1 worker slice

    # load spatial data for objects in the FOV
    masterLoadStage = lsst.ap.pipeline.LoadStage(0, loadPolicy)
    masterLoadStage.initialize(mq1, mq0)
    masterLoadStage.setRun(runId)
    masterLoadStage.setRank(1)
    masterLoadStage.setUniverseSize(2)
    workerLoadStage = lsst.ap.pipeline.LoadStage(0, loadPolicy)
    workerLoadStage.initialize(wq1, wq0)
    workerLoadStage.setRun(runId)
    workerLoadStage.setRank(0)
    workerLoadStage.setUniverseSize(2)

    # read difference sources, match, write match results
    masterMatchStage1Input = lsst.dps.IOStage.InputStage(1, match1InputPolicy)
    masterMatchStage1Input.initialize(mq2, mq1)
    masterMatchStage1Input.setRun(runId)
    masterMatchStage1Input.setRank(1)
    masterMatchStage1Input.setUniverseSize(2)
    workerMatchStage1Input = lsst.dps.IOStage.InputStage(1, match1InputPolicy)
    workerMatchStage1Input.initialize(wq2, wq1)
    workerMatchStage1Input.setRun(runId)
    workerMatchStage1Input.setRank(0)
    workerMatchStage1Input.setUniverseSize(2)

    masterMatchStage1 = lsst.ap.pipeline.MatchDiaSourcesStage(2, None)
    masterMatchStage1.initialize(mq3, mq2)
    masterMatchStage1.setRun(runId)
    masterMatchStage1.setRank(1)
    masterMatchStage1.setUniverseSize(2)
    workerMatchStage1 = lsst.ap.pipeline.MatchDiaSourcesStage(2, None)
    workerMatchStage1.initialize(wq3, wq2)
    workerMatchStage1.setRun(runId)
    workerMatchStage1.setRank(0)
    workerMatchStage1.setUniverseSize(2)

    masterMatchStage1Output = lsst.dps.IOStage.OutputStage(3, match1OutputPolicy)
    masterMatchStage1Output.initialize(mq4, mq3)
    masterMatchStage1Output.setRun(runId)
    masterMatchStage1Output.setRank(1)
    masterMatchStage1Output.setUniverseSize(2)
    workerMatchStage1Output = lsst.dps.IOStage.OutputStage(3, match1OutputPolicy)
    workerMatchStage1Output.initialize(wq4, wq3)
    workerMatchStage1Output.setRun(runId)
    workerMatchStage1Output.setRank(0)
    workerMatchStage1Output.setUniverseSize(2)

    # read moving object predictions, match, write match results
    masterMatchStage2Input = lsst.dps.IOStage.InputStage(4, match2InputPolicy)
    masterMatchStage2Input.initialize(mq5, mq4)
    masterMatchStage2Input.setRun(runId)
    masterMatchStage2Input.setRank(1)
    masterMatchStage2Input.setUniverseSize(2)
    workerMatchStage2Input = lsst.dps.IOStage.InputStage(4, match2InputPolicy)
    workerMatchStage2Input.initialize(wq5, wq4)
    workerMatchStage2Input.setRun(runId)
    workerMatchStage2Input.setRank(0)
    workerMatchStage2Input.setUniverseSize(2)

    masterMatchStage2 = lsst.ap.pipeline.MatchMopsPredsStage(5, None)
    masterMatchStage2.initialize(mq6, mq5)
    masterMatchStage2.setRun(runId)
    masterMatchStage2.setRank(1)
    masterMatchStage2.setUniverseSize(2)
    workerMatchStage2 = lsst.ap.pipeline.MatchMopsPredsStage(5, None)
    workerMatchStage2.initialize(wq6, wq5)
    workerMatchStage2.setRun(runId)
    workerMatchStage2.setRank(0)
    workerMatchStage2.setUniverseSize(2)

    masterMatchStage2Output = lsst.dps.IOStage.OutputStage(6, match2OutputPolicy)
    masterMatchStage2Output.initialize(mq7, mq6)
    masterMatchStage2Output.setRun(runId)
    masterMatchStage2Output.setRank(1)
    masterMatchStage2Output.setUniverseSize(2)
    workerMatchStage2Output = lsst.dps.IOStage.OutputStage(6, match2OutputPolicy)
    workerMatchStage2Output.initialize(wq7, wq6)
    workerMatchStage2Output.setRun(runId)
    workerMatchStage2Output.setRank(0)
    workerMatchStage2Output.setUniverseSize(2)

    # store spatial data for new objects in the FOV
    masterStoreStage = lsst.ap.pipeline.StoreStage(7, storePolicy)
    masterStoreStage.initialize(mq8, mq7)
    masterStoreStage.setRun(runId)
    masterStoreStage.setRank(1)
    masterStoreStage.setUniverseSize(2)
    workerStoreStage = lsst.ap.pipeline.StoreStage(7, storePolicy)
    workerStoreStage.initialize(wq8, wq7)
    workerStoreStage.setRun(runId)
    workerStoreStage.setRank(0)
    workerStoreStage.setUniverseSize(2)

    # Create the TCS event triggering a single visit
    dt        = lsst.mwi.persistence.DateTime(long(time.time())*1000000000)
    visitTime = dt.utc2mjd()
    triggerAssociationEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitId', 1))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitTime', visitTime))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('filterName', 'u'))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVRA', 273.48066298343))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVDec', -27.125))

    # Create the event triggering the match against moving object predictions
    triggerMatchMopsPredsEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    triggerMatchMopsPredsEvent.addProperty(lsst.mwi.data.DataProperty('visitId', 1))


    # Create clipboards for the master and worker, add the events to both,
    # and place them onto the input queues
    masterClipboard = lsst.dps.Clipboard.Clipboard()
    workerClipboard = lsst.dps.Clipboard.Clipboard()

    masterClipboard.put('triggerAssociationEvent', triggerAssociationEvent)
    workerClipboard.put('triggerAssociationEvent', triggerAssociationEvent)

    masterClipboard.put('triggerMatchMopsPredsEvent', triggerMatchMopsPredsEvent)
    workerClipboard.put('triggerMatchMopsPredsEvent', triggerMatchMopsPredsEvent)

    mq0.addDataset(masterClipboard)
    wq0.addDataset(workerClipboard)

    # Now run it!
    try:
        # load
        masterLoadStage.preprocess()
        workerLoadStage.process()
        masterLoadStage.postprocess()
        
        # match difference sources
        masterMatchStage1Input.preprocess()
        workerMatchStage1Input.process()
        masterMatchStage1Input.postprocess()
        
        masterMatchStage1.preprocess()
        workerMatchStage1.process()
        masterMatchStage1.postprocess()
        
        masterMatchStage1Output.preprocess()
        workerMatchStage1Output.process()
        masterMatchStage1Output.postprocess()
        
        # match moving object predictions
        masterMatchStage2Input.preprocess()
        workerMatchStage2Input.process()
        masterMatchStage2Input.postprocess()
        
        masterMatchStage2.preprocess()
        workerMatchStage2.process()
        masterMatchStage2.postprocess()
        
        masterMatchStage2Output.preprocess()
        workerMatchStage2Output.process()
        masterMatchStage2Output.postprocess()
        
        # store
        masterStoreStage.preprocess()
        workerStoreStage.process()
        masterStoreStage.postprocess()
    except Exception, e:
        print e

    # Close log to void bogus memory-leak reports
    lsst.mwi.logging.Log.closeDefaultLog()
    TestCleanup.cleanup()


if __name__ == '__main__':
    runOneVisit()
    # check for memory leaks
    if lsst.mwi.data.Citizen_census(0, 0) != 0:
        print lsst.mwi.data.Citizen_census(0, 0), 'Objects leaked:'
        print lsst.mwi.data.Citizen_census(lsst.mwi.data.cout, 0)

