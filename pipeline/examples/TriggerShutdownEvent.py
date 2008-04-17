#! /usr/bin/env python

import lsst.daf.base
import lsst.ctrl.events

if __name__ == "__main__":
    print 'starting...\n'
    externalEventTransmitter = lsst.ctrl.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'shutdownAssociationEvent')
    shutdownAssociationEvent = lsst.daf.base.DataProperty.createPropertyNode('root')
    externalEventTransmitter.publish('eventtype', shutdownAssociationEvent)

