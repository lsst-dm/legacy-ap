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

#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/coord/Coord.h"

#include "../Common.h"


namespace lsst { namespace ap { namespace utils {

/** Clamps the given latitude/declination to <tt>[-M_PI/2, M_PI/2]</tt>.
  */
inline lsst::afw::geom::Angle const clampPhi(lsst::afw::geom::Angle const a) {
    using lsst::afw::geom::Angle;
    using lsst::afw::geom::radians;
    using lsst::afw::geom::HALFPI;

    if (static_cast<double>(a) < -HALFPI) {
        return Angle(-HALFPI, radians);
    } else if (static_cast<double>(a) > HALFPI) {
        return Angle(HALFPI, radians);
    }
    return a;
}

void thetaRangeReduce(lsst::afw::geom::Angle &min, lsst::afw::geom::Angle &max);

lsst::afw::geom::Angle const maxAlpha(lsst::afw::geom::Angle radius,
                                      lsst::afw::geom::Angle centerPhi);

void positionAndVelocity(Eigen::Vector3d &p,
                         Eigen::Vector3d &v,
                         lsst::afw::geom::Angle ra,
                         lsst::afw::geom::Angle decl,
                         double muRa,
                         double muDecl,
                         double vRadial,
                         lsst::afw::geom::Angle parallax);

Eigen::Vector2d const cartesianToSpherical(Eigen::Vector3d const &v);

inline lsst::afw::coord::IcrsCoord const cartesianToIcrs(Eigen::Vector3d const &v) {
    Eigen::Vector2d sc = cartesianToSpherical(v);
    return lsst::afw::coord::IcrsCoord(sc.x() * lsst::afw::geom::radians,
                                       sc.y() * lsst::afw::geom::radians);
}

/** Returns the angular separation (rad) between two 3-vectors of arbitrary
  * magnitude.
  */
inline lsst::afw::geom::Angle const angularSeparation(
    Eigen::Vector3d const &v1,
    Eigen::Vector3d const &v2)
{
    using lsst::afw::geom::Angle;
    using lsst::afw::geom::radians;

    double ss = v1.cross(v2).norm();
    double cs = v1.dot(v2);
    if (ss != 0.0 || cs != 0.0) {
        return Angle(std::atan2(ss, cs), radians);
    }
    return Angle(0.0, radians);
}

}}} // namespace lsst::ap::utils

#endif // LSST_AP_UTILS_SPATIALUTILS_H

