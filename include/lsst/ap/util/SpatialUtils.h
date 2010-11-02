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
  * @brief  Spatial utility functions.
  */
#ifndef LSST_AP_UTIL_SPATIALUTILS_H
#define LSST_AP_UTIL_SPATIALUTILS_H

#include "lsst/afw/geom/Point.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace util {

/** Converts degrees to radians.
  */
inline double degrees(double const rad) {
    return rad*DEGREES_PER_RADIAN;
}

/** Converts radians to degrees.
  */
inline double radians(double const deg) {
    return deg*RADIANS_PER_DEGREE;
}

/** Clamps the given latitude/declination to [-90, 90].
  */
inline double clampPhi(double const phi) {
    return phi <= -90.0 ? -90.0 : (phi >= 90.0 ? 90.0 : phi);
}

void thetaRangeReduce(double &min, double &max);

double maxAlpha(double radius, double centerPhi);

lsst::afw::geom::Point3D const starProperMotion(double ra,
                                                double decl,
                                                double muRa,
                                                double muDecl,
                                                double vRadial,
                                                double parallax,
                                                double fromEpoch,
                                                double toEpoch);

lsst::afw::geom::Point2D const cartesianToSpherical(lsst::afw::geom::Point3D const &v);

}}} // namespace lsst::ap::util

#endif // LSST_AP_UTIL_SPATIALUTILS_H
