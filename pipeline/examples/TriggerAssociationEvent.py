#! /usr/bin/env python

import lsst.mwi.data
import lsst.mwi.persistence
import lsst.events
import time

if __name__ == "__main__":
    print 'starting...\n'
    externalEventTransmitter = lsst.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerAssociationEvent')

    dt        = lsst.mwi.persistence.DateTime(long(time.time())*1000000000)
    visitTime = dt.utc2mjd()
    triggerAssociationEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitId', 708125))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitTime', visitTime))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('filterName', 'u'))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVRA', 273.48066298343))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVDec', -27.125))

    externalEventTransmitter.publish('eventtype', triggerAssociationEvent)

