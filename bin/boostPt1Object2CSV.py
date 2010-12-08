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

import csv
import math
import optparse
from textwrap import dedent

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

def rangeReduce(val, nullstr):
    if val == None:
        return nullstr
    v = math.fmod(math.degrees(val), 360.0)
    if v < 0.0:
        v += 360.0
        if v == 360.0:
            v = 0.0
    return v

def deg(val, nullstr):
    if val == None or val == nullstr:
        return nullstr
    return math.degrees(val)

def deg2(val, nullstr):
    if val == None or val == nullstr:
        return nullstr
    return math.degrees(math.degrees(val))


def convertFilter(object, row, filter, nullstr):
    i = filterMap[filter]
    if object.hasFilter(i):
        pfa = object.getPerFilterAttributes(i)
        row.append(pfa.getNumObs())
        row.extend([nullstr]*14)
        row.append(pfa.getFlux() or nullstr)
        row.append(pfa.getFluxSigma() or nullstr)
        row.extend([nullstr]*5)
        row.append(pfa.getEarliestObsTime())
        row.append(pfa.getLatestObsTime())
        row.extend([nullstr, nullstr])
        row.append(pfa.getE1() or nullstr)
        row.append(pfa.getE1Sigma() or nullstr)
        row.append(pfa.getE2() or nullstr)
        row.append(pfa.getE2Sigma() or nullstr)
        row.append(pfa.getRadius() or nullstr)
        row.append(pfa.getRadiusSigma() or nullstr)
        row.append(pfa.getFlags())
    else:
        row.extend([nullstr]*33)

def objects2CSV(objects, csvWriter, nullstr):
    for o in objects:
        row = []
        row.append(o.getClusterId())
        row.append(nullstr)
        row.append(rangeReduce(o.getRaPs(), nullstr))
        row.append(deg(o.getRaPsSigma(), nullstr))
        row.append(deg(o.getPsDec(), nullstr))
        row.append(deg(o.getPsDecSigma(), nullstr))
        row.append(deg2(o.getRaDecPsCov(), nullstr))
        row.append(rangeReduce(o.getRaSg(), nullstr))
        row.append(deg(o.getRaSgSigma(), nullstr))
        row.append(deg(o.getSgDec(), nullstr))
        row.append(deg(o.getSgDecSigma(), nullstr))
        row.append(deg2(o.getRaDecSgCov(), nullstr))
        row.extend([nullstr]*12)
        row.append(o.getEarliestObsTime())
        row.append(o.getLatestObsTime())
        row.append(o.getFlags())
        for filter in ('u', 'g', 'r', 'i', 'z', 'y'):
            convertFilter(o, row, filter, nullstr)
        csvWriter.writerow(row)


if __name__ == "__main__":
    usage = """\
usage: %prog [options] in_file_1 in_file_2 ... out_csv_file

    Reads in one or more lsst.ap.cluster.PersistableSourceClusterVector
    objects persisted using BoostStorage. These are converted
    and stored to a CSV file suitable for qserv ingest.
    """
    parser = optparse.OptionParser(usage)
    # CSV format options
    fmt = optparse.OptionGroup(
        parser, "CSV format options", dedent("""\
        See http://docs.python.org/library/csv.html#csv-fmt-params for
        details."""))
    fmt.add_option(
        "-D", "--delimiter", dest="delimiter", default=",",
        help=dedent("""\
        One character string used to separate fields in the
        input CSV files. The default is %default."""))
    fmt.add_option(
        "-n", "--no-doublequote", dest="doublequote", action="store_false",
        help=dedent("""\
        Turn off double quoting of quote characters inside a CSV field."""))
    fmt.add_option(
        "-e", "--escapechar", dest="escapechar", default=None,
        help="Delimiter escape character.")
    quoteHelp = dedent("""\
        CSV quoting style. May be one of %d  (quote all fields), %d (quote
        fields containing special characters), %d (quote non-numeric fields)
        or %d (never quote fields). The default is %%default.""" %
        (csv.QUOTE_ALL, csv.QUOTE_MINIMAL,
         csv.QUOTE_NONNUMERIC, csv.QUOTE_NONE))
    fmt.add_option(
        "-Q", "--quoting", type="int", dest="quoting",
        default=csv.QUOTE_MINIMAL, help=quoteHelp)
    fmt.add_option(
        "-q", "--quotechar", dest="quotechar", default='"',
        help=dedent("""\
        One character string to quote fields with; defaults to %default."""))
    fmt.add_option(
        "-N", "--null", dest="null", default=r"\N", help=dedent("""\
        String to output for null fields."""))
    parser.add_option_group(fmt)
    (opts, inputs) = parser.parse_args()
    if len(inputs) < 2:
        parser.error("Please specify at least one input file and an output file")
    outFile = open(inputs.pop(), 'wb+')
    outWriter = csv.writer(outFile,
                           delimiter=opts.delimiter,
                           doublequote=opts.doublequote,
                           escapechar=opts.escapechar,
                           quoting=opts.quoting,
                           quotechar=opts.quotechar,
                           lineterminator='\n')
    persistence = dafPersistence.Persistence.getPersistence(pexPolicy.Policy())
    props = dafBase.PropertySet()
    for f in inputs:
        sl = dafPersistence.StorageList()
        loc = dafPersistence.LogicalLocation(f)
        sl.append(persistence.getRetrieveStorage("BoostStorage", loc))
        persistable = persistence.unsafeRetrieve(
            "PersistableSourceClusterVector", sl, props)
        psv = apCluster.PersistableSourceClusterVector.swigConvert(persistable)
        objects2CSV(psv.getClusters(), outWriter, opts.null)
 
