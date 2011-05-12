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
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_UTILS_SPATIALUTILS_H
#define LSST_AP_UTILS_SPATIALUTILS_H

#include <math.h>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "../Common.h"


namespace lsst { namespace ap { namespace utils {

/** Clamps the given latitude/declination to <tt>[-M_PI/2, M_PI/2]</tt>.
  */
inline double clampPhi(double const a) {
    return a <= -M_PI_2 ? -M_PI_2 : (a >= M_PI_2 ? M_PI_2 : a);
}

LSST_AP_API void thetaRangeReduce(double &min, double &max);

LSST_AP_API double maxAlpha(double radius, double centerPhi);

LSST_AP_API void positionAndVelocity(Eigen::Vector3d &p,
                                     Eigen::Vector3d &v,
                                     double ra,
                                     double decl,
                                     double muRa,
                                     double muDecl,
                                     double vRadial,
                                     double parallax);

LSST_AP_API Eigen::Vector2d const cartesianToSpherical(Eigen::Vector3d const &v);

/** Converts spherical coordinate <tt>(theta, phi)</tt> (rad) to a unit 3-vector.
  */
inline Eigen::Vector3d const sphericalToCartesian(double ra, double dec) {
    double cosDec = std::cos(dec);
    return Eigen::Vector3d(std::cos(ra)*cosDec,
                           std::sin(ra)*cosDec,
                           std::sin(dec));
}

/** Returns the angular separation (rad) between two 3-vectors of arbitrary
  * magnitude.
  */
inline double angularSeparation(Eigen::Vector3d const &v1,
                                Eigen::Vector3d const &v2) {
    double ss = v1.cross(v2).norm();
    double cs = v1.dot(v2);
    return (ss != 0.0 || cs != 0.0) ? std::atan2(ss, cs) : 0.0;
}

}}} // namespace lsst::ap::utils

#endif // LSST_AP_UTILS_SPATIALUTILS_H

