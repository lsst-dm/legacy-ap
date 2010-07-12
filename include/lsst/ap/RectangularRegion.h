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
 * @brief   Class describing a rectangular (in ra and dec) region on the sky.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_RECTANGULAR_REGION_H
#define LSST_AP_RECTANGULAR_REGION_H

#include "Common.h"
#include "CircularRegion.h"


namespace lsst { namespace ap {

/** @brief  A rectangular region (in right ascension and declination) of the unit sphere. */
class LSST_AP_API RectangularRegion {

public :

    RectangularRegion() : _minRa(0.0), _maxRa(0.0), _minDec(0.0), _maxDec(0.0) {}

    RectangularRegion(
        double const minRa,
        double const maxRa,
        double const minDec,
        double const maxDec
    );

    RectangularRegion(
        double const cenRa,
        double const cenDec,
        double const radius
    );

    explicit RectangularRegion(CircularRegion const & region);

    /// Returns the minimum right ascension of points in the region.
    double getMinRa() const {
        return _minRa;
    }

    /// Returns the maximum right ascension of points in the region.
    double getMaxRa() const {
        return _maxRa;
    }

    /// Returns the minimum declination of points in the region
    double getMinDec() const {
        return _minDec;
    }

    /// Returns the maximum declination of points in the region
    double getMaxDec() const {
        return _maxDec;
    }

private :

    double _minRa;
    double _maxRa;
    double _minDec;
    double _maxDec;

    void fromCircle(double const centerRa, double const centerDec, double const radius);
};


}}  // end of namespace lsst::ap

#endif // LSST_AP_RECTANGULAR_REGION_H
