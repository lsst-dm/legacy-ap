#! /usr/bin/env python

import lsst.mwi.data
import lsst.events
import datetime

if __name__ == "__main__":
    print 'starting...\n'
    externalEventTransmitter = lsst.events.EventTransmitter('lsst8.ncsa.uiuc.edu', 'triggerAssociationEvent')

    visitTime = datetime.datetime.utcnow().isoformat(' ')
    triggerAssociationEvent = lsst.mwi.data.SupportFactory.createPropertyNode('root')
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitId', 1))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('visitTime', visitTime))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('filterName', 'u'))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVRA', 273.48066298343))
    triggerAssociationEvent.addProperty(lsst.mwi.data.DataProperty('FOVDec', -27.125))

    externalEventTransmitter.publish('eventtype', triggerAssociationEvent)

