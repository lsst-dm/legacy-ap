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


def translate(source, method, flag=-1):
    val = method(source)
    if source.isNull(flag) or val != val:
        return None
    return val

def translateDeg(source, method, flag=-1):
    val = translate(source, method, flag)
    if val != None:
        val = math.degrees(val)
    return val

def translateDegRangeReduce(source, method, flag=-1):
    val = translate(source, method, flag)
    if val != None:
        val = math.fmod(math.degrees(val), 360.0)
        if val < 0.0:
            val += 360.0
    return val

def sources2CSV(sources, csvWriter):
    for s in sources:
        row = []
        row.append(translate(s, Source.getSourceId))
        row.append(translate(s, Source.getAmpExposureId, afwDetection.AMP_EXPOSURE_ID))
        row.append(translate(s, Source.getFilterId))
        row.append(translate(s, Source.getObjectId, afwDetection.OBJECT_ID))
        row.append(translate(s, Source.getMovingObjectId, afwDetection.MOVING_OBJECT_ID))
        row.append(translate(s, Source.getProcHistoryId))
        row.append(translateDegRangeReduce(s, Source.getRa))
        row.append(translateDeg(s, Source.getRaErrForDetection, afwDetection.RA_ERR_FOR_DETECTION))
        row.append(translateDeg(s, Source.getRaErrForWcs))
        row.append(translateDeg(s, Source.getDec))
        row.append(translateDeg(s, Source.getDecErrForDetection, afwDetection.DEC_ERR_FOR_DETECTION))
        row.append(translateDeg(s, Source.getDecErrForWcs))
        row.append(translate(s, Source.getXFlux, afwDetection.X_FLUX))
        row.append(translate(s, Source.getXFluxErr, afwDetection.X_FLUX_ERR))
        row.append(translate(s, Source.getYFlux, afwDetection.Y_FLUX))
        row.append(translate(s, Source.getYFluxErr, afwDetection.Y_FLUX_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaFlux, afwDetection.RA_FLUX))
        row.append(translateDeg(s, Source.getRaFluxErr, afwDetection.RA_FLUX_ERR))
        row.append(translateDeg(s, Source.getDecFlux, afwDetection.DEC_FLUX))
        row.append(translateDeg(s, Source.getDecFluxErr, afwDetection.DEC_FLUX_ERR))
        row.append(translate(s, Source.getXPeak, afwDetection.X_PEAK))
        row.append(translate(s, Source.getYPeak, afwDetection.Y_PEAK))
        row.append(translateDegRangeReduce(s, Source.getRaPeak, afwDetection.RA_PEAK))
        row.append(translateDeg(s, Source.getDecPeak, afwDetection.DEC_PEAK))
        row.append(translate(s, Source.getXAstrom, afwDetection.X_ASTROM))
        row.append(translate(s, Source.getXAstromErr, afwDetection.X_ASTROM_ERR))
        row.append(translate(s, Source.getYAstrom, afwDetection.Y_ASTROM))
        row.append(translate(s, Source.getYAstromErr, afwDetection.Y_ASTROM_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaAstrom, afwDetection.RA_ASTROM))
        row.append(translateDeg(s, Source.getRaAstromErr, afwDetection.RA_ASTROM_ERR))
        row.append(translateDeg(s, Source.getDecAstrom, afwDetection.DEC_ASTROM))
        row.append(translateDeg(s, Source.getDecAstromErr, afwDetection.DEC_ASTROM_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaObject, afwDetection.RA_OBJECT))
        row.append(translateDeg(s, Source.getDecObject, afwDetection.DEC_OBJECT))
        row.append(translate(s, Source.getTaiMidPoint))
        row.append(translate(s, Source.getTaiRange, afwDetection.TAI_RANGE))
        row.append(translate(s, Source.getPsfFlux))
        row.append(translate(s, Source.getPsfFluxErr))
        row.append(translate(s, Source.getApFlux))
        row.append(translate(s, Source.getApFluxErr))
        row.append(translate(s, Source.getModelFlux))
        row.append(translate(s, Source.getModelFluxErr))
        row.append(translate(s, Source.getXFlux, afwDetection.X_FLUX))
        row.append(translate(s, Source.getXFluxErr, afwDetection.X_FLUX_ERR))
        row.append(translate(s, Source.getYFlux, afwDetection.Y_FLUX))
        row.append(translate(s, Source.getYFluxErr, afwDetection.Y_FLUX_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaFlux, afwDetection.RA_FLUX))
        row.append(translateDeg(s, Source.getRaFluxErr, afwDetection.RA_FLUX_ERR))
        row.append(translateDeg(s, Source.getDecFlux, afwDetection.DEC_FLUX))
        row.append(translateDeg(s, Source.getDecFluxErr, afwDetection.DEC_FLUX_ERR))
        row.append(translate(s, Source.getXPeak, afwDetection.X_PEAK))
        row.append(translate(s, Source.getYPeak, afwDetection.Y_PEAK))
        row.append(translateDegRangeReduce(s, Source.getRaPeak, afwDetection.RA_PEAK))
        row.append(translateDeg(s, Source.getDecPeak, afwDetection.DEC_PEAK))
        row.append(translate(s, Source.getXAstrom, afwDetection.X_ASTROM))
        row.append(translate(s, Source.getXAstromErr, afwDetection.X_ASTROM_ERR))
        row.append(translate(s, Source.getYAstrom, afwDetection.Y_ASTROM))
        row.append(translate(s, Source.getYAstromErr, afwDetection.Y_ASTROM_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaAstrom, afwDetection.RA_ASTROM))
        row.append(translateDeg(s, Source.getRaAstromErr, afwDetection.RA_ASTROM_ERR))
        row.append(translateDeg(s, Source.getDecAstrom, afwDetection.DEC_ASTROM))
        row.append(translateDeg(s, Source.getDecAstromErr, afwDetection.DEC_ASTROM_ERR))
        row.append(translateDegRangeReduce(s, Source.getRaObject, afwDetection.RA_OBJECT))
        row.append(translateDeg(s, Source.getDecObject, afwDetection.DEC_OBJECT))
        row.append(translate(s, Source.getTaiMidPoint))
        row.append(translate(s, Source.getTaiRange, afwDetection.TAI_RANGE))
        row.append(translate(s, Source.getPsfFlux))
        row.append(translate(s, Source.getPsfFluxErr))
        row.append(translate(s, Source.getApFlux))
        row.append(translate(s, Source.getApFluxErr))
        row.append(translate(s, Source.getModelFlux))
        row.append(translate(s, Source.getPetroFlux, afwDetection.PETRO_FLUX))
        row.append(translate(s, Source.getPetroFluxErr, afwDetection.PETRO_FLUX_ERR))
        row.append(translate(s, Source.getInstFlux))
        row.append(translate(s, Source.getInstFluxErr))
        row.append(translate(s, Source.getNonGrayCorrFlux, afwDetection.NON_GRAY_CORR_FLUX))
        row.append(translate(s, Source.getNonGrayCorrFluxErr, afwDetection.NON_GRAY_CORR_FLUX_ERR))
        row.append(translate(s, Source.getAtmCorrFlux, afwDetection.ATM_CORR_FLUX))
        row.append(translate(s, Source.getAtmCorrFluxErr, afwDetection.ATM_CORR_FLUX_ERR))
        row.append(translate(s, Source.getApDia, afwDetection.AP_DIA))
        row.append(translate(s, Source.getIxx, afwDetection.IXX))
        row.append(translate(s, Source.getIxxErr, afwDetection.IXX_ERR))
        row.append(translate(s, Source.getIyy, afwDetection.IYY))
        row.append(translate(s, Source.getIyyErr, afwDetection.IYY_ERR))
        row.append(translate(s, Source.getIxy, afwDetection.IXY))
        row.append(translate(s, Source.getIxyErr, afwDetection.IXY_ERR))
        row.append(translate(s, Source.getSnr))
        row.append(translate(s, Source.getChi2))
        row.append(translate(s, Source.getSky, afwDetection.SKY))
        row.append(translate(s, Source.getSkyErr, afwDetection.SKY_ERR))
        row.append(translate(s, Source.getFlagForAssociation, afwDetection.FLAG_FOR_ASSOCIATION))
        row.append(translate(s, Source.getFlagForDetection, afwDetection.FLAG_FOR_DETECTION))
        row.append(translate(s, Source.getFlagForWcs))
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
        persistable = persistence.unsafeRetrieve("PersistableSourceVector", sl, props)
        psv = afwDetection.PersistableSourceVector.swigConvert(persistable)
        sources2CSV(psv.getSources(), outWriter)
 
