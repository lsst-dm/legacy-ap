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
#include "lsst/ap/match/ImageInfo.h"

#include "boost/make_shared.hpp"

#include "lsst/daf/base/DateTime.h"


using lsst::daf::base::DateTime;


namespace lsst { namespace ap { namespace match {

ImageInfo::ImageInfo(
    int64_t id, ///< Image identifier, typically a science CCD ID
    lsst::daf::base::PropertySet::Ptr metadata ///< FITS metadata for image
) :
    _id(id),
    _extent(),
    _calib(boost::make_shared<lsst::afw::image::Calib>()),
    _wcs(lsst::afw::image::makeWcs(metadata, false)),
    _canCalibrateFlux(false)
{
    // get image extents
    _extent.setX(metadata->getAsInt("NAXIS1"));
    _extent.setY(metadata->getAsInt("NAXIS2"));
    // get image exposure time and mid-point
    DateTime midPoint(metadata->getAsDouble("TIME-MID"),
                      DateTime::MJD, DateTime::TAI);
    _calib->setMidTime(midPoint);
    _calib->setExptime(metadata->getAsDouble("EXPTIME"));
    // get image flux
    if (metadata->exists("FLUXMAG0") && metadata->exists("FLUXMAG0ERR")) {
        _calib->setFluxMag0(metadata->getAsDouble("FLUXMAG0"),
                            metadata->getAsDouble("FLUXMAG0ERR"));
        _canCalibrateFlux = true;
    }
}

ImageInfo::~ImageInfo() { }

}}} // namespace lsst::ap::match

