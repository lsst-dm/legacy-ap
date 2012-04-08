// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

/** @file
  * @brief  ImageInfo implementation
  * @author Serge Monkewitz
  */
#include "lsst/ap/match/ExposureInfo.h"

#include <algorithm>
#include <string>

#include "boost/algorithm/string/trim.hpp"
#include "boost/regex.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/daf/base/DateTime.h"
#include "lsst/afw/image/Filter.h"
#include "lsst/afw/image/ImageUtils.h"

#include "lsst/ap/utils/Csv.h"
#include "lsst/ap/utils/SpatialUtils.h"

using std::max;

using boost::algorithm::trim_copy;

using lsst::pex::exceptions::InvalidParameterException;
using lsst::pex::exceptions::RuntimeErrorException;
using lsst::pex::logging::Log;
using lsst::daf::base::DateTime;
using lsst::daf::base::PropertySet;
using lsst::afw::image::Filter;
using lsst::afw::image::PixelZeroPos;

using lsst::ap::utils::angularSeparation;
using lsst::ap::utils::cartesianToIcrs;
using lsst::ap::utils::CsvControl;
using lsst::ap::utils::CsvReader;
using lsst::ap::utils::maxAlpha;


namespace lsst { namespace ap { namespace match {

// -- ExposureInfo implementation ----

std::string const ExposureInfo::DEF_ID_KEY("scienceCcdExposureId");

ExposureInfo::ExposureInfo(
    lsst::daf::base::PropertySet::Ptr metadata, ///< FITS metadata for image
    std::string const &idKey ///< Metadata key for unique integer image id
) :
    _center(),
    _earthPos(),
    _id(metadata->getAsInt64(idKey)),
    _fluxMag0(std::numeric_limits<double>::quiet_NaN()),
    _fluxMag0Sigma(std::numeric_limits<double>::quiet_NaN()),
    _extent(),
    _wcs(lsst::afw::image::makeWcs(metadata, false)),
    _filter(trim_copy(metadata->getAsString("FILTER")), false),
    _canCalibrateFlux(false),
    _epValid(false)
{
    // get image extents
    _extent.setX(metadata->getAsInt("NAXIS1"));
    _extent.setY(metadata->getAsInt("NAXIS2"));
    // get image exposure time and mid-point
    DateTime mid(metadata->getAsString("TIME-MID"));
    _epoch = mid.get(DateTime::MJD, DateTime::TAI);
    _expTime = metadata->getAsDouble("EXPTIME");
    // get image flux
    if (metadata->exists("FLUXMAG0")) {
        if (metadata->exists("FLUXMAG0ERR")) {
            _fluxMag0Sigma = metadata->getAsDouble("FLUXMAG0ERR");
            if (_fluxMag0Sigma < 0.0) {
                throw LSST_EXCEPT(InvalidParameterException,
                                  "negative magnitude zero point error");
            }
        } else {
            _fluxMag0Sigma = 0.0;
        }
        _fluxMag0 = metadata->getAsDouble("FLUXMAG0");
        if (_fluxMag0 <= 0.0) {
            throw LSST_EXCEPT(InvalidParameterException,
                              "magnitude zero point is negative or zero");
        }
        _canCalibrateFlux = true;
    }
    // compute image center and corners
    Eigen::Vector3d c = _pixToSky(0.5*_extent.getX() + PixelZeroPos,
                                  0.5*_extent.getY() + PixelZeroPos);
    Eigen::Vector3d llc = _pixToSky(PixelZeroPos - 0.5,
                                    PixelZeroPos - 0.5);
    Eigen::Vector3d ulc = _pixToSky(PixelZeroPos - 0.5,
                                    _extent.getY() - 0.5 + PixelZeroPos);
    Eigen::Vector3d lrc = _pixToSky(_extent.getX() - 0.5 + PixelZeroPos,
                                    PixelZeroPos - 0.5);
    Eigen::Vector3d urc = _pixToSky(_extent.getX() - 0.5 + PixelZeroPos,
                                    _extent.getY() - 0.5 + PixelZeroPos);
    // compute bounding box from bounding circle
    _center = cartesianToIcrs(c);
    _radius = angularSeparation(c, llc);
    _radius = max(_radius, angularSeparation(c, ulc));
    _radius = max(_radius, angularSeparation(c, lrc));
    _radius = max(_radius, angularSeparation(c, urc));
    _alpha = maxAlpha(_radius, _center.getLatitude());
}

ExposureInfo::~ExposureInfo() { }

/** Returns a flux value calibrated using the flux of a zero magnitude object
  * associated with this exposure.
  */
double ExposureInfo::calibrateFlux(
    double flux,            ///< flux to calibrate, DN
    double fluxScale        ///< flux scaling factor, must be \> 0.0
) const {
    if (fluxScale <= 0.0) {
        throw LSST_EXCEPT(InvalidParameterException, "flux "
                          "scaling factor is <= 0");
    }
    if (!canCalibrateFlux()) {
        throw LSST_EXCEPT(RuntimeErrorException, "Cannot calibrate "
                          "flux without fluxMag0");
    }
    return fluxScale*flux/_fluxMag0;
}

/** Returns a calibrated flux and its variance. The flux of a zero magnitude
  * object associated with this exposure is used to perform the calibration.
  */
std::pair<double, double> const ExposureInfo::calibrateFlux(
    double flux,            ///< flux to calibrate, DN
    double fluxSigma,       ///< standard deviation of @a flux
    double fluxScale        ///< flux scaling factor, must be \> 0.0
) const {
    if (fluxScale <= 0.0) {
        throw LSST_EXCEPT(InvalidParameterException, "flux "
                          "scaling factor is <= 0");
    }
    if (fluxSigma <= 0.0) {
        throw LSST_EXCEPT(InvalidParameterException, "flux error is <= 0");
    }
    if (!canCalibrateFlux()) {
        throw LSST_EXCEPT(RuntimeErrorException, "Cannot calibrate "
                          "flux without fluxMag0");
    }
    // Use delta method to estimate the variance oF flux/fluxMag0.
    double f02 = _fluxMag0*_fluxMag0;
    double f0Var = _fluxMag0Sigma*_fluxMag0Sigma;
    double f2 = flux*flux;
    double fVar = fluxSigma*fluxSigma;
    double var = (fVar*f02 + f0Var*f2)/(f02*f02);
    return std::make_pair(fluxScale*flux/_fluxMag0, fluxScale*fluxScale*var);
}

double ExposureInfo::getMinCoord0() const {
    return static_cast<double>(_center.getLongitude() - _alpha);
}

double ExposureInfo::getMaxCoord0() const {
    return static_cast<double>(_center.getLongitude() + _alpha);
}

double ExposureInfo::getMinCoord1() const {
    return static_cast<double>(_center.getLatitude() - _radius);
}

double ExposureInfo::getMaxCoord1() const {
    return static_cast<double>(_center.getLatitude() + _radius);
}

Eigen::Vector3d const ExposureInfo::_pixToSky(double x, double y) const {
    return _wcs->pixelToSky(x, y)->toIcrs().getVector().asEigen();
}


// -- ExposureInfoMap implementation ----

ExposureInfoMap::ExposureInfoMap() : _map() { }

ExposureInfoMap::~ExposureInfoMap() { }

void ExposureInfoMap::insert(ExposureInfo::Ptr info) {
    if (!info) {
        throw LSST_EXCEPT(InvalidParameterException, "Cannot insert a "
                          "null ExposureInfo into an ExposureInfoMap");
    }
    int64_t id = info->getId();
    if (contains(id)) {
        throw LSST_EXCEPT(InvalidParameterException, "ExposureInfoMap "
                          "already contains an exposure with the "
                          "specified id");
    }
    _map.insert(std::pair<int64_t, ExposureInfo::Ptr>(id, info));
}

void ExposureInfoMap::clear() {
    _map.clear();
}

bool ExposureInfoMap::erase(int64_t id) {
    return _map.erase(id) != 0u;
}


namespace {

int const INT_T = 0;
int const DOUBLE_T = 1;
int const STRING_T = 2;

std::tr1::unordered_map<std::string, int> const & fitsKeyMap() {
    static char const * const I_KEYS[] = {
        "NAXIS1", "NAXIS2",
        "A_ORDER", "B_ORDER",
        "AP_ORDER", "BP_ORDER"
    };
    static char const * const D_KEYS[] = {
        "EQUINOX", "EXPTIME",
        "FLUXMAG0", "FLUXMAG0ERR",
        "CRPIX1", "CRPIX2",
        "CRVAL1", "CRVAL2",
        "CDELT1", "CDELT2",
        "CD1_1", "CD1_2", "CD2_1", "CD2_2",
        "PC1_1", "PC1_2", "PC2_1", "PC2_2"
    };
    static char const * const S_KEYS[] = {
        "RADESYS", "TIME-MID", "FILTER",
        "CUNIT1", "CUNIT2",
        "CTYPE1", "CTYPE2"
    };
    static std::tr1::unordered_map<std::string, int> map;
    if (map.size() == 0) {
        // build key to type mapping for the FITS cards we care about
        for (size_t i = 0; i < sizeof(I_KEYS)/sizeof(char const *); ++i) {
            map.insert(std::make_pair(std::string(I_KEYS[i]), INT_T));
        }
        for (size_t i = 0; i < sizeof(D_KEYS)/sizeof(char const *); ++i) {
            map.insert(std::make_pair(std::string(D_KEYS[i]), DOUBLE_T));
        }
        for (size_t i = 0; i < sizeof(S_KEYS)/sizeof(char const *); ++i) {
            map.insert(std::make_pair(std::string(S_KEYS[i]), STRING_T));
        }
    }
    return map;
}

} // namespace lsst::ap::match::<anonymous>


/** Reads an exposure metadata key-value CSV file (where metadata keys
  * must have been grouped by exposure id). An ExposureInfo object
  * is created for each input exposure and appended to @a exposures.
  */
void readExposureInfos(
    std::vector<ExposureInfo::Ptr> & exposures, ///< ExposureInfo objects are
                                                ///  appended to this vector.
    std::string const & csvFile,                ///< Metadata table path.
    CsvControl  const & control,                ///< Metadata table CSV format.
    std::string const & idColumn                ///< Name of ID column, e.g. "scienceCcdExposureId".
) {
    typedef std::tr1::unordered_map<std::string, int> FkMap;
    typedef FkMap::const_iterator FkIter;

    static boost::regex const sipRegex("^[AB]P?_[0-9]+_[0-9]+$");

    Log log(Log::getDefaultLog(), "lsst.ap.match");
    log.log(Log::INFO, "Reading exposure metadata from " + csvFile);

    FkMap const & fkMap = fitsKeyMap();
    // open CSV file and get column indexes
    CsvReader reader(csvFile, control, true);
    int const idCol = reader.getIndexOf(idColumn);
    int const keyCol = reader.getIndexOf("metadataKey");
    int const intCol = reader.getIndexOf("intValue");
    int const doubleCol = reader.getIndexOf("doubleValue");
    int const stringCol = reader.getIndexOf("stringValue");
    if (idCol == -1) {
        throw LSST_EXCEPT(RuntimeErrorException, "Exposure metadata table "
                          "has no column named " + idCol);
    }
    if (keyCol == -1) {
        throw LSST_EXCEPT(RuntimeErrorException, "Exposure metadata table "
                          "has no column named metadataKey");
    }
    if (intCol == -1) {
        throw LSST_EXCEPT(RuntimeErrorException, "Exposure metadata table "
                          "has no column named intValue");
    }
    if (doubleCol == -1) {
        throw LSST_EXCEPT(RuntimeErrorException, "Exposure metadata table "
                          "has no column named doubleValue");
    }
    if (stringCol == -1) {
        throw LSST_EXCEPT(RuntimeErrorException, "Exposure metadata table "
                          "has no column named stringValue");
    }
    PropertySet::Ptr ps;
    int64_t lastId = reader.get<int64_t>(idCol) - 1;
    for (; !reader.isDone(); reader.nextRecord()) {
        int64_t id = reader.get<int64_t>(idCol);
        if (id != lastId) {
            if (ps) {
                exposures.push_back(ExposureInfo::Ptr(new ExposureInfo(ps)));
            }
            ps.reset(new PropertySet());
            ps->set(idColumn, id);
            lastId = id;
        }
        std::string const key = reader.get(keyCol);
        FkIter k = fkMap.find(key);
        if (k != fkMap.end()) {
            if (k->second == INT_T) { // integer-valued key
                ps->set(key, reader.get<int>(intCol));
            } else if (k->second == DOUBLE_T) { // double-valued key
                ps->set(key, reader.get<double>(doubleCol));
            } else { // string-valued key
                ps->set(key, reader.get(stringCol));
            }
        } else if (boost::regex_search(key.begin(), key.end(), sipRegex)) {
            // key is a SIP coefficient
            ps->set(key, reader.get<double>(doubleCol));
        }
    }
    if (ps) {
        exposures.push_back(ExposureInfo::Ptr(new ExposureInfo(ps)));
    }
}

}}} // namespace lsst::ap::match

