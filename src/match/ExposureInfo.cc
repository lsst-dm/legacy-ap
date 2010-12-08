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
#include "boost/make_shared.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/daf/base/DateTime.h"
#include "lsst/afw/image/ImageUtils.h"

#include "lsst/ap/utils/SpatialUtils.h"

using std::max;

using boost::algorithm::trim_copy;

using lsst::pex::exceptions::InvalidParameterException;
using lsst::pex::exceptions::RuntimeErrorException;
using lsst::daf::base::DateTime;
using lsst::afw::image::PixelZeroPos;

using lsst::ap::utils::cartesianToSpherical;
using lsst::ap::utils::angularSeparation;
using lsst::ap::utils::maxAlpha;


namespace lsst { namespace ap { namespace match {

// -- ExposureInfo implementation ----

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
    _canCalibrateFlux(false),
    _epValid(false)
{
    std::string filter = trim_copy(metadata->getAsString("FILTER"));
    if (filter == "u") {
        _filterId = 0;
    } else if (filter == "g") {
        _filterId = 1;
    } else if (filter == "r") {
        _filterId = 2;
    } else if (filter == "i") {
        _filterId = 3;
    } else if (filter == "z") {
        _filterId = 4;
    } else if (filter == "y" || filter == "i2") {
        _filterId = 5;
    } else {
        throw LSST_EXCEPT(InvalidParameterException,
                          "unknown filter " + filter);
    }
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
    _center = cartesianToSpherical(c);
    _radius = angularSeparation(c, llc);
    _radius = max(_radius, angularSeparation(c, ulc));
    _radius = max(_radius, angularSeparation(c, lrc));
    _radius = max(_radius, angularSeparation(c, urc));
    _alpha = maxAlpha(_radius, _center.y());
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
    return _center.x() - _alpha;
}

double ExposureInfo::getMaxCoord0() const {
    return _center.x() + _alpha;
}

double ExposureInfo::getMinCoord1() const {
    return _center.y() - _radius;
}

double ExposureInfo::getMaxCoord1() const {
    return _center.y() + _radius;
}

Eigen::Vector3d const ExposureInfo::_pixToSky(double x, double y) const {
    return _wcs->pixelToSky(x, y)->toIcrs().getVector().asVector();
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


}}} // namespace lsst::ap::match

