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
 * @brief   Class describing a circular region on the sky.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_CIRCULAR_REGION_H
#define LSST_AP_CIRCULAR_REGION_H

#include <algorithm>

#include "Common.h"


namespace lsst { namespace ap {

/** @brief  A circular region of the unit sphere (sky). */
class LSST_AP_API CircularRegion {

public :

    CircularRegion() : _centerRa(0.0), _centerDec(0.0), _radius(0.0) {}

    CircularRegion(
        double const ra,
        double const dec,
        double const radius
    );

    /// Returns the right ascension of the circle center.
    double getCenterRa() const {
        return _centerRa;
    }

    /// Returns the declination of the circle center.
    double getCenterDec() const {
        return _centerDec;
    }

    /// Returns the radius of the circle.
    double getRadius() const {
        return _radius;
    }

    /// Returns the minimum declination of points in the region
    double getMinDec() const {
        double d = _centerDec - _radius;
        return (d <= -90.0 ? -90.0 : d);
    }

    /// Returns the maximum declination of points in the region
    double getMaxDec() const {
        double d = _centerDec + _radius;
        return (d >= 90.0 ? 90.0 : d);
    }

private :

    double _centerRa;
    double _centerDec;
    double _radius;
};


}}  // end of namespace lsst::ap

#endif // LSST_AP_CIRCULAR_REGION_H
