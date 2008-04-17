#! /usr/bin/env python

import lsst.daf.base
import lsst.ctrl.events
import time

if __name__ == "__main__":
    print 'starting...\n'
    evTransmitter1 = lsst.ctrl.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerAssociationEvent')
    evTransmitter2 = lsst.ctrl.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerMatchMopsPredsEvent')
    evTransmitter3 = lsst.ctrl.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'shutdownAssociationEvent')

    dt        = lsst.daf.base.DateTime(long(time.time())*1000000000)
    visitTime = dt.utc2mjd()
    triggerAssociationEvent = lsst.daf.base.DataProperty.createPropertyNode('root')
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('visitId', 708125))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('visitTime', visitTime))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('filterName', 'u'))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('FOVRA', 333.880166667))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('FOVDec', -17.7374166667))

    triggerMatchMopsPredsEvent = lsst.daf.base.DataProperty.createPropertyNode('root')
    triggerMatchMopsPredsEvent.addProperty(lsst.daf.base.DataProperty('visitId', 708125))

    for i in xrange(1000):
        print "Starting iteration %i" % i
        evTransmitter1.publish('eventtype', triggerAssociationEvent)
        time.sleep(10)
        evTransmitter2.publish('eventtype', triggerMatchMopsPredsEvent)
        time.sleep(10)

    shutdownAssociationEvent = lsst.daf.base.DataProperty.createPropertyNode('root')
    evTransmitter3.publish('eventtype', shutdownAssociationEvent)
