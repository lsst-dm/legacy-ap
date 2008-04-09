#! /usr/bin/env python

import lsst.daf.base
import lsst.ctrl.events

if __name__ == "__main__":
    print 'starting...\n'
    externalEventTransmitter = lsst.ctrl.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerMatchMopsPredsEvent')

    triggerMatchMopsPredsEvent = lsst.daf.base.DataProperty.createPropertyNode('root')
    triggerMatchMopsPredsEvent.addProperty(lsst.daf.base.DataProperty('visitId', 708125))

    externalEventTransmitter.publish('eventtype', triggerMatchMopsPredsEvent)

