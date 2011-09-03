#! /usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

from __future__ import with_statement
import math
import optparse

import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.daf.persistence as dafPersistence
import lsst.ap.cluster as apCluster

filterMap = { 'u': 0,
              'g': 1,
              'r': 2,
              'i': 3,
              'z': 4,
              'y': 5, 'i2': 5
            }

nullstr = '\\N'
inf = float('INF')

def isSpecial(val):
    return val == None or val != val or val == inf or val == -inf

def rangeReduce(val):
    if isSpecial(val):
        return nullstr
    v = math.fmod(math.degrees(val), 360.0)
    if v < 0.0:
        v += 360.0
        if v == 360.0:
            v = 0.0
    return repr(v)

def deg(val):
    if isSpecial(val):
        return nullstr
    return repr(math.degrees(val))

def arcsec(val):
    if isSpecial(val):
        return nullstr
    return repr(math.degrees(val)*3600.0)

def deg2(val):
    if isSpecial(val):
        return nullstr
    return repr(math.degrees(math.degrees(val)))

def translate(val):
    if isSpecial(val):
        return nullstr
    return repr(val)

def convertFilter(object, row, filter):
    i = filterMap[filter]
    if object.hasFilter(i):
        pfa = object.getPerFilterAttributes(i)
        row.append(translate(pfa.getNumObs()))
        row.extend([nullstr]*14)
        row.append(translate(pfa.getPsFlux()))
        row.append(translate(pfa.getPsFluxSigma()))
        row.append(translate(pfa.getSgFlux()))
        row.append(translate(pfa.getSgFluxSigma()))
        row.append(translate(pfa.getGaussianFlux()))
        row.append(translate(pfa.getGaussianFluxSigma()))
        row.append(nullstr)
        row.append(translate(pfa.getEarliestObsTime()))
        row.append(translate(pfa.getLatestObsTime()))
        row.extend([nullstr, nullstr])
        row.append(translate(pfa.getE1()))
        row.append(translate(pfa.getE1Sigma()))
        row.append(translate(pfa.getE2()))
        row.append(translate(pfa.getE2Sigma()))
        row.append(arcsec(pfa.getRadius()))
        row.append(arcsec(pfa.getRadiusSigma()))
        row.append(translate(pfa.getNumPsFluxSamples()))
        row.append(translate(pfa.getNumSgFluxSamples()))
        row.append(translate(pfa.getNumGaussianFluxSamples()))
        row.append(translate(pfa.getNumEllipticitySamples()))
        row.append(nullstr)
    else:
        row.extend([nullstr]*37)

def objects2CSV(objects, csvWriter):
    for o in objects:
        row = []
        row.append(translate(o.getClusterId()))
        row.append(nullstr)
        row.append(rangeReduce(o.getRaPs()))
        row.append(deg(o.getRaPsSigma()))
        row.append(deg(o.getDecPs()))
        row.append(deg(o.getDecPsSigma()))
        row.append(deg2(o.getRaDecPsCov()))
        row.append(rangeReduce(o.getRaSg()))
        row.append(deg(o.getRaSgSigma()))
        row.append(deg(o.getDecSg()))
        row.append(deg(o.getDecSgSigma()))
        row.append(deg2(o.getRaDecSgCov()))
        row.extend([nullstr]*12)
        row.append(translate(o.getEarliestObsTime()))
        row.append(translate(o.getLatestObsTime()))
        row.append(translate(o.getMeanObsTime()))
        row.append(translate(o.getFlags()))
        for filter in ('u', 'g', 'r', 'i', 'z', 'y'):
            convertFilter(o, row, filter)
        csvWriter.write(','.join(row))
        csvWriter.write('\n')


if __name__ == "__main__":
    usage = """\
usage: %prog [options] in_file_1 in_file_2 ... out_csv_file

    Reads in one or more lsst.ap.cluster.PersistableSourceClusterVector
    objects persisted using BoostStorage. These are converted
    and stored to a CSV file suitable for qserv ingest.
    """
    parser = optparse.OptionParser(usage)
    (opts, inputs) = parser.parse_args()
    if len(inputs) < 2:
        parser.error("Please specify at least one input file and an output file")
    persistence = dafPersistence.Persistence.getPersistence(pexPolicy.Policy())
    props = dafBase.PropertySet()
    with open(inputs.pop(), 'wb+') as csvWriter:
        for f in inputs:
            sl = dafPersistence.StorageList()
            loc = dafPersistence.LogicalLocation(f)
            sl.append(persistence.getRetrieveStorage("BoostStorage", loc))
            persistable = persistence.unsafeRetrieve(
                "PersistableSourceClusterVector", sl, props)
            psv = apCluster.PersistableSourceClusterVector.swigConvert(persistable)
            objects2CSV(psv.getClusters(), csvWriter)
 
