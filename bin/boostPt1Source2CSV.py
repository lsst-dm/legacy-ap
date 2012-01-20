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
import lsst.afw.detection as afwDetection
# without this, unpersisting sources can fail because of
# unregistered classes. In the future, a meas_multifit
# import will become necessary as well.
import lsst.meas.algorithms as measAlgorithms

from lsst.afw.detection import Source

nullstr = '\\N'
inf = float('INF')

def isSpecial(val):
    return val == None or val != val or val == inf or val == -inf

def translate(source, method, flag=-1):
    val = method(source)
    if source.isNull(flag) or isSpecial(val):
        return nullstr
    return repr(val)

def translateDeg(source, method, flag=-1):
    val = method(source)
    val = val.asDegrees()
    if source.isNull(flag) or isSpecial(val):
        return nullstr
    return repr(val)

def translateDegRangeReduce(source, method, flag=-1):
    val = method(source)
    val = math.fmod(val.asDegrees(), 360.0)
    if source.isNull(flag) or isSpecial(val):
        return nullstr
    if val < 0.0:
        val += 360.0
        if val == 360.0:
            val = 0.0
    return repr(val)

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
        row.append(translate(s, Source.getPsfIxx, afwDetection.PSF_IXX))
        row.append(translate(s, Source.getPsfIxxErr, afwDetection.PSF_IXX_ERR))
        row.append(translate(s, Source.getPsfIyy, afwDetection.PSF_IYY))
        row.append(translate(s, Source.getPsfIyyErr, afwDetection.PSF_IYY_ERR))
        row.append(translate(s, Source.getPsfIxy, afwDetection.PSF_IXY))
        row.append(translate(s, Source.getPsfIxyErr, afwDetection.PSF_IXY_ERR))
        row.append(translate(s, Source.getE1, afwDetection.E1))
        row.append(translate(s, Source.getE1Err, afwDetection.E1_ERR))
        row.append(translate(s, Source.getE2, afwDetection.E2))
        row.append(translate(s, Source.getE2Err, afwDetection.E2_ERR))
        row.append(translate(s, Source.getResolution, afwDetection.RESOLUTION))
        row.append(translate(s, Source.getShear1, afwDetection.SHEAR1))
        row.append(translate(s, Source.getShear1Err, afwDetection.SHEAR1_ERR))
        row.append(translate(s, Source.getShear2, afwDetection.SHEAR2))
        row.append(translate(s, Source.getShear2Err, afwDetection.SHEAR2_ERR))
        row.append(translate(s, Source.getSigma, afwDetection.SIGMA))
        row.append(translate(s, Source.getSigmaErr, afwDetection.SIGMA_ERR))
        row.append(translate(s, Source.getShapeStatus, afwDetection.SHAPE_STATUS))
        row.append(translate(s, Source.getSnr))
        row.append(translate(s, Source.getChi2))
        row.append(translate(s, Source.getSky, afwDetection.SKY))
        row.append(translate(s, Source.getSkyErr, afwDetection.SKY_ERR))
        row.append(translate(s, Source.getFlagForAssociation))
        row.append(translate(s, Source.getFlagForDetection))
        row.append(translate(s, Source.getFlagForWcs))
        csvWriter.write(','.join(row))
        csvWriter.write('\n')

if __name__ == "__main__":
    usage = """\
usage: %prog [options] in_file_1 in_file_2 ... out_csv_file

    Reads in one or more lsst.afw.detection.PersistableSourceVectors
    persisted using BoostStorage from the given files. These are converted
    and stored to a CSV file suitable for qserv ingest.
    """
    parser = optparse.OptionParser(usage)
    (opts, inputs) = parser.parse_args()
    if len(inputs) < 2:
        parser.error("Please specify at least one input file and an output file")
    persistence = dafPersistence.Persistence.getPersistence(pexPolicy.Policy())
    props = dafBase.PropertySet()
    with open(inputs.pop(), 'wb') as csvWriter:
        for f in inputs:
            sl = dafPersistence.StorageList()
            loc = dafPersistence.LogicalLocation(f)
            sl.append(persistence.getRetrieveStorage("BoostStorage", loc))
            persistable = persistence.unsafeRetrieve(
                "PersistableSourceVector", sl, props)
            psv = afwDetection.PersistableSourceVector.swigConvert(persistable)
            sources2CSV(psv.getSources(), csvWriter)
 
