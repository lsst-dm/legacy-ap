#! /usr/bin/env python
import csv
import math
import optparse
from textwrap import dedent

import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.daf.persistence as dafPersistence
import lsst.afw.detection as afwDetection

from lsst.afw.detection import Source


def translate(source, method, nullstr, flag=-1):
    val = method(source)
    if source.isNull(flag) or val != val:
        return nullstr
    return val

def translateDeg(source, method, nullstr, flag=-1):
    val = translate(source, method, flag)
    if val != nullstr:
        val = math.degrees(val)
    return val

def translateDegRangeReduce(source, method, nullstr, flag=-1):
    val = translate(source, method, flag)
    if val != nullstr:
        val = math.fmod(math.degrees(val), 360.0)
        if val < 0.0:
            val += 360.0
    return val

def sources2CSV(sources, csvWriter, nullstr):
    for s in sources:
        row = []
        row.append(translate(s, Source.getSourceId, nullstr))
        row.append(translate(s, Source.getAmpExposureId, nullstr, afwDetection.AMP_EXPOSURE_ID))
        row.append(translate(s, Source.getFilterId, nullstr))
        row.append(translate(s, Source.getObjectId, nullstr, afwDetection.OBJECT_ID))
        row.append(translate(s, Source.getMovingObjectId, nullstr, afwDetection.MOVING_OBJECT_ID))
        row.append(translate(s, Source.getProcHistoryId, nullstr))
        row.append(translateDegRangeReduce(s, Source.getRa, nullstr))
        row.append(translateDeg(s, Source.getRaErrForDetection, nullstr, afwDetection.RA_ERR_FOR_DETECTION))
        row.append(translateDeg(s, Source.getRaErrForWcs, nullstr))
        row.append(translateDeg(s, Source.getDec, nullstr))
        row.append(translateDeg(s, Source.getDecErrForDetection, nullstr, afwDetection.DEC_ERR_FOR_DETECTION))
        row.append(translateDeg(s, Source.getDecErrForWcs, nullstr))
        row.append(translate(s, Source.getXFlux, nullstr, afwDetection.X_FLUX))
        row.append(translate(s, Source.getXFluxErr, nullstr, afwDetection.X_FLUX_ERR))
        row.append(translate(s, Source.getYFlux, nullstr, afwDetection.Y_FLUX))
        row.append(translate(s, Source.getYFluxErr, nullstr, afwDetection.Y_FLUX_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaFlux, nullstr, afwDetection.RA_FLUX))
        row.append(translateDeg(s, Source.getRaFluxErr, nullstr, afwDetection.RA_FLUX_ERR))
        row.append(translateDeg(s, Source.getDecFlux, nullstr, afwDetection.DEC_FLUX))
        row.append(translateDeg(s, Source.getDecFluxErr, nullstr, afwDetection.DEC_FLUX_ERR))
        row.append(translate(s, Source.getXPeak, nullstr, afwDetection.X_PEAK))
        row.append(translate(s, Source.getYPeak, nullstr, afwDetection.Y_PEAK))
        row.append(translateDegRangeReduce(s, Source.getRaPeak, nullstr, afwDetection.RA_PEAK))
        row.append(translateDeg(s, Source.getDecPeak, nullstr, afwDetection.DEC_PEAK))
        row.append(translate(s, Source.getXAstrom, nullstr, afwDetection.X_ASTROM))
        row.append(translate(s, Source.getXAstromErr, nullstr, afwDetection.X_ASTROM_ERR))
        row.append(translate(s, Source.getYAstrom, nullstr, afwDetection.Y_ASTROM))
        row.append(translate(s, Source.getYAstromErr, nullstr, afwDetection.Y_ASTROM_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaAstrom, nullstr, afwDetection.RA_ASTROM))
        row.append(translateDeg(s, Source.getRaAstromErr, nullstr, afwDetection.RA_ASTROM_ERR))
        row.append(translateDeg(s, Source.getDecAstrom, nullstr, afwDetection.DEC_ASTROM))
        row.append(translateDeg(s, Source.getDecAstromErr, nullstr, afwDetection.DEC_ASTROM_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaObject, nullstr, afwDetection.RA_OBJECT))
        row.append(translateDeg(s, Source.getDecObject, nullstr, afwDetection.DEC_OBJECT))
        row.append(translate(s, Source.getTaiMidPoint, nullstr))
        row.append(translate(s, Source.getTaiRange, nullstr, afwDetection.TAI_RANGE))
        row.append(translate(s, Source.getPsfFlux, nullstr))
        row.append(translate(s, Source.getPsfFluxErr, nullstr))
        row.append(translate(s, Source.getApFlux, nullstr))
        row.append(translate(s, Source.getApFluxErr, nullstr))
        row.append(translate(s, Source.getModelFlux, nullstr))
        row.append(translate(s, Source.getModelFluxErr, nullstr))
        row.append(translate(s, Source.getPetroFlux, nullstr, afwDetection.PETRO_FLUX))
        row.append(translate(s, Source.getPetroFluxErr, nullstr, afwDetection.PETRO_FLUX_ERR))
        row.append(translate(s, Source.getInstFlux, nullstr))
        row.append(translate(s, Source.getInstFluxErr, nullstr))
        row.append(translate(s, Source.getNonGrayCorrFlux, nullstr, afwDetection.NON_GRAY_CORR_FLUX))
        row.append(translate(s, Source.getNonGrayCorrFluxErr, nullstr, afwDetection.NON_GRAY_CORR_FLUX_ERR))
        row.append(translate(s, Source.getAtmCorrFlux, nullstr, afwDetection.ATM_CORR_FLUX))
        row.append(translate(s, Source.getAtmCorrFluxErr, nullstr, afwDetection.ATM_CORR_FLUX_ERR))
        row.append(translate(s, Source.getApDia, nullstr, afwDetection.AP_DIA))
        row.append(translate(s, Source.getIxx, nullstr, afwDetection.IXX))
        row.append(translate(s, Source.getIxxErr, nullstr, afwDetection.IXX_ERR))
        row.append(translate(s, Source.getIyy, nullstr, afwDetection.IYY))
        row.append(translate(s, Source.getIyyErr, nullstr, afwDetection.IYY_ERR))
        row.append(translate(s, Source.getIxy, nullstr, afwDetection.IXY))
        row.append(translate(s, Source.getIxyErr, nullstr, afwDetection.IXY_ERR))
        row.append(translate(s, Source.getSnr, nullstr))
        row.append(translate(s, Source.getChi2, nullstr))
        row.append(translate(s, Source.getSky, nullstr, afwDetection.SKY))
        row.append(translate(s, Source.getSkyErr, nullstr, afwDetection.SKY_ERR))
        row.append(translate(s, Source.getFlagForAssociation, nullstr, afwDetection.FLAG_FOR_ASSOCIATION))
        row.append(translate(s, Source.getFlagForDetection, nullstr, afwDetection.FLAG_FOR_DETECTION))
        row.append(translate(s, Source.getFlagForWcs, nullstr))
        csvWriter.writerow(row)

if __name__ == "__main__":
    usage = """\
usage: %prog [options] in_file_1 in_file_2 ... out_csv_file

    Reads in one or more lsst.afw.detection.PersistableSourceVectors
    persisted using BoostStorage from the given files. These are converted
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
            "PersistableSourceVector", sl, props)
        psv = afwDetection.PersistableSourceVector.swigConvert(persistable)
        sources2CSV(psv.getSources(), outWriter, opts.null)
 
