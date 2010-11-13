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

#include "boost/make_shared.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/daf/base/DateTime.h"

#include "lsst/ap/utils/SpatialUtils.h"

using std::max;

using lsst::pex::exceptions::InvalidParameterException;
using lsst::daf::base::DateTime;
using lsst::ap::utils::cartesianToSpherical;
using lsst::ap::utils::angularSeparation;
using lsst::ap::utils::maxAlpha;


namespace lsst { namespace ap { namespace match {

ExposureInfo::ExposureInfo(
    lsst::daf::base::PropertySet::Ptr metadata, ///< FITS metadata for image
    std::string const &idKey ///< Metadata key for unique integer image id
) :
    _center(),
    _earthPos(),
    _id(metadata->getAsInt64(idKey)),
    _fluxMag0(std::numeric_limits<double>::quiet_NaN()),
    _fluxMag0Err(std::numeric_limits<double>::quiet_NaN()),
    _extent(),
    _wcs(lsst::afw::image::makeWcs(metadata, false)),
    _canCalibrateFlux(false),
    _epValid(false)
{
    std::string filter = metadata->getAsString("FILTER");
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
    _epoch = metadata->getAsDouble("TIME-MID");
    _expTime = metadata->getAsDouble("EXPTIME");
    // get image flux
    if (metadata->exists("FLUXMAG0") && metadata->exists("FLUXMAG0ERR")) {
        _fluxMag0 = metadata->getAsDouble("FLUXMAG0");
        _fluxMag0Err = metadata->getAsDouble("FLUXMAG0ERR");
        _canCalibrateFlux = true;
    }
    // compute image center and corners
    Eigen::Vector3d c = _pixToSky(0.5*_extent.getX(), 0.5*_extent.getY());
    Eigen::Vector3d llc = _pixToSky(-0.5, -0.5);
    Eigen::Vector3d ulc = _pixToSky(-0.5, _extent.getY() - 0.5);
    Eigen::Vector3d lrc = _pixToSky(_extent.getX() - 0.5, -0.5);
    Eigen::Vector3d urc = _pixToSky(_extent.getX() - 0.5, _extent.getY() - 0.5);
    // compute bounding box from bounding circle
    _center = cartesianToSpherical(c);
    _radius = angularSeparation(c, llc);
    _radius = max(_radius, angularSeparation(c, ulc));
    _radius = max(_radius, angularSeparation(c, lrc));
    _radius = max(_radius, angularSeparation(c, urc));
    _alpha = maxAlpha(_radius, _center.y());
}

ExposureInfo::~ExposureInfo() { }

Eigen::Vector3d const ExposureInfo::_pixToSky(double x, double y) const {
    return _wcs->pixelToSky(x, y)->toIcrs().getVector().asVector();
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

}}} // namespace lsst::ap::match

