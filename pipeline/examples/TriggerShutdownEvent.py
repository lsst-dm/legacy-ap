#! /usr/bin/env python

import lsst.mwi.data
import lsst.mwi.persistence
import lsst.events

if __name__ == "__main__":
    print 'starting...\n'
    externalEventTransmitter = lsst.events.EventTransmitter('lsst4.ncsa.uiuc.edu', 'shutdownAssociationEvent')
    shutdownAssociationEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    externalEventTransmitter.publish('eventtype', shutdownAssociationEvent)

