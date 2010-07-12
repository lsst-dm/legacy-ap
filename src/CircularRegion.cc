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
 

/**
 * @file
 * @brief   Implementation of the CircularRegion class.
 *
 * @ingroup ap
 */

#include "lsst/pex/exceptions.h"

#include "lsst/ap/CircularRegion.h"


namespace ex = lsst::pex::exceptions;

lsst::ap::CircularRegion::CircularRegion(
    double const ra,
    double const dec,
    double const radius
) :
    _centerRa(ra),
    _centerDec(dec),
    _radius(radius)
{
    if (ra < 0.0 || ra >= 360.0) {
        throw LSST_EXCEPT(ex::RangeErrorException,
                          "right ascension must be in range [0, 360) degrees");
    }
    if (dec < -90.0 || dec > 90.0) {
        throw LSST_EXCEPT(ex::RangeErrorException,
                          "declination must be in range  [-90, 90] degrees");
    }
    if (radius < 0.0 || radius > 90.0) {
        throw LSST_EXCEPT(ex::RangeErrorException,
                          "circle radius must be in range  [0, 90] degrees");
    }
}

