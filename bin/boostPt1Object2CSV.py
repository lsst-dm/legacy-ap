#! /usr/bin/env python
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

def rangeReduce(val):
    if val == None:
        return None
    v = math.fmod(math.degrees(val), 360.0)
    if v < 0.0:
        v += 360.0
    return v

def deg(val):
    if val == None:
        return None
    return math.degrees(val)


def convertFilter(object, row, filter):
    i = filterMap[filter]
    if object.hasFilter(i):
        pfa = object.getPerFilterAttributes(i)
        row.append(pfa.getNumObs())
        row.extend([None]*14)
        row.append(pfa.getFlux())
        row.append(pfa.getFluxSigma())
        row.extend([None]*5)
        row.append(pfa.getEarliestObsTime())
        row.append(pfa.getLatestObsTime())
        row.extend([None, None])
        row.append(pfa.getE1())
        row.append(pfa.getE1Sigma())
        row.append(pfa.getE2())
        row.append(pfa.getE2Sigma())
        row.append(pfa.getRadius())
        row.append(pfa.getRadiusSigma())
        row.append(pfa.getFlags())
    else:
        row.extend([None]*33)

def objects2CSV(objects, csvWriter):
    for o in objects:
        row = []
        row.append(o.getClusterId())
        row.append(None)
        row.append(rangeReduce(o.getRa()))
        row.append(deg(o.getRaSigma()))
        row.append(deg(o.getDec()))
        row.append(deg(o.getDecSigma()))
        row.append(deg(deg(o.getRaDecCov())))
        row.extend([None]*17)
        row.append(o.getEarliestObsTime())
        row.append(o.getLatestObsTime())
        row.append(o.getFlags())
        for filter in ('u', 'g', 'r', 'i', 'z', 'y'):
            convertFilter(o, row, filter)
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
        "-I", "--skipinitialspace", dest="skipinitialspace",
        action="store_true", help=dedent("""\
        Ignore whitespace immediately following delimiters."""))
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
                           skipinitialspace=opts.skipinitialspace,
                           lineterminator='\n')
    persistence = dafPersistence.Persistence.getPersistence(pexPolicy.Policy())
    props = dafBase.PropertySet()
    for f in inputs:
        sl = dafPersistence.StorageList()
        loc = dafPersistence.LogicalLocation(f)
        sl.append(persistence.getRetrieveStorage("BoostStorage", loc))
        persistable = persistence.unsafeRetrieve("PersistableSourceClusterVector", sl, props)
        psv = apCluster.PersistableSourceClusterVector.swigConvert(persistable)
        objects2CSV(psv.getClusters(), outWriter)
 
