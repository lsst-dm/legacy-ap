#! /usr/bin/env python

import lsst.mwi.data
import lsst.mwi.persistence
import lsst.events
import time

if __name__ == "__main__":
    print 'starting...\n'
    evTransmitter1 = lsst.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerAssociationEvent')
    evTransmitter2 = lsst.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerMatchMopsPredsEvent')
    evTransmitter3 = lsst.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'shutdownAssociationEvent')

    dt        = lsst.mwi.persistence.DateTime(long(time.time())*1000000000)
    visitTime = dt.utc2mjd()
    triggerAssociationEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitId', 708125))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitTime', visitTime))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('filterName', 'u'))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVRA', 333.880166667))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVDec', -17.7374166667))

    triggerMatchMopsPredsEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    triggerMatchMopsPredsEvent.addProperty(lsst.mwi.data.DataProperty('visitId', 708125))

    for i in xrange(1000):
        print "Starting iteration %i" % i
        evTransmitter1.publish('eventtype', triggerAssociationEvent)
        time.sleep(10)
        evTransmitter2.publish('eventtype', triggerMatchMopsPredsEvent)
        time.sleep(10)

    shutdownAssociationEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    evTransmitter3.publish('eventtype', shutdownAssociationEvent)
