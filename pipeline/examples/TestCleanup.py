#! /usr/bin/env python
"""
Cleans up after a run of AssociationPipelineTest.py or an MPI based run of the association pipeline.

Run with:
   python TestCleanup.py
"""

import lsst.daf.persistence


def swallow(function):
    """Calls the given function and swallows any exception it might throw"""
    try:
        function()
    except:
        pass

def cleanup():
    """Cleans up after a test run of the association pipeline"""
    loc = lsst.daf.persistence.LogicalLocation('mysql://lsst10.ncsa.uiuc.edu:3306/test')
    db  = lsst.daf.persistence.DbStorage()
    db.setPersistLocation(loc)
    db.startTransaction()
    #swallow(lambda : db.dropTable('DiaSourceToObjectMatches_visit708125'))
    #swallow(lambda : db.dropTable('MopsPredToDiaSourceMatches_visit708125'))
    swallow(lambda : db.dropTable('NewObjectIdPairs_visit708125'))
    swallow(lambda : db.dropTable('BestMatch_visit708125'))
    swallow(lambda : db.dropTable('MatchedObjects_visit708125'))
    swallow(lambda : db.dropTable('NewObjects_visit708125'))
    swallow(lambda : db.truncateTable('DIASource'))
    swallow(lambda : db.truncateTable('VarObject'))
    db.endTransaction()

if __name__ == "__main__":
    cleanup()

