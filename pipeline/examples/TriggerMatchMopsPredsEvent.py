#! /usr/bin/env python

import lsst.mwi.data
import lsst.events

if __name__ == "__main__":
    print 'starting...\n'
    externalEventTransmitter = lsst.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerMatchMopsPredsEvent')

    triggerMatchMopsPredsEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    triggerMatchMopsPredsEvent.addProperty(lsst.mwi.data.DataProperty('visitId', 1))

    externalEventTransmitter.publish('eventtype', triggerMatchMopsPredsEvent)

