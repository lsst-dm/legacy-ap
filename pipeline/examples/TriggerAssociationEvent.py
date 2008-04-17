#! /usr/bin/env python

import lsst.daf.base
import lsst.ctrl.events
import time

if __name__ == "__main__":
    print 'starting...\n'
    externalEventTransmitter = lsst.ctrl.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerAssociationEvent')

    dt        = lsst.daf.base.DateTime(long(time.time())*1000000000)
    visitTime = dt.utc2mjd()
    triggerAssociationEvent = lsst.daf.base.DataProperty.createPropertyNode('root')
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('visitId', 708125))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('visitTime', visitTime))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('filterName', 'u'))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('FOVRA', 333.880166667))
    triggerAssociationEvent.addProperty(lsst.daf.base.DataProperty('FOVDec', -17.7374166667))

    externalEventTransmitter.publish('eventtype', triggerAssociationEvent)

