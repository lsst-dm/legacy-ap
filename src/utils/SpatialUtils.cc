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
  * @brief  Implementation of spatial utility functions.
  */
#include "lsst/ap/utils/SpatialUtils.h"

#include <cfloat>

#include "lsst/pex/exceptions.h"
#include "lsst/ap/constants.h"

using lsst::pex::exceptions::InvalidParameterError;
using lsst::pex::exceptions::RuntimeError;

using lsst::afw::geom::HALFPI;
using lsst::afw::geom::PI;
using lsst::afw::geom::TWOPI;
using lsst::afw::geom::Angle;
using lsst::afw::geom::radians;
using lsst::afw::coord::IcrsCoord;


namespace lsst { namespace ap { namespace utils {

/** Reduces a theta (longitude/right-ascension) range. The resulting
  * range will have <tt> min \> max </tt> if it wraps across the 0/2*PI
  * radiandiscontinuity. Valid inputs are:
  *
  * @li  any <tt> min, max </tt> with <tt> min \<= max </tt>
  * @li  <tt> min \> max </tt> so long as
  *      <tt> min \<= 2*PI && max \>= 0.0 </tt>
  */
void thetaRangeReduce(Angle &min, Angle &max) {
    if (min > max) {
        if (max < 0.0 || min >= TWOPI) {
            throw LSST_EXCEPT(InvalidParameterError,
                              "Invalid longitude angle interval");
        }
    } else if (max - min >= TWOPI) {
        min = 0.0 * radians;
        max = TWOPI * radians;
    } else if (min < 0.0 || max >= TWOPI) {
        // range reduce
        min.wrap();
        max.wrap();
    }
}

/** Computes the extent in longitude angle <tt> [-alpha, alpha] </tt> of
  * the circle with radius @a radius and center at latitude angle @a centerPhi.
  *
  * @pre  <tt> radius >= 0.0 && radius <= PI/2 </tt>
  *
  * Note that @a centerPhi is clamped to lie in [-PI/2, PI/2].
  */
Angle const maxAlpha(
    Angle radius,   ///< circle radius
    Angle centerPhi ///< latitude angle of circle center
) {
    static const double POLE_EPSILON = 1e-7;

    if (radius < 0.0 || radius > HALFPI) {
        throw LSST_EXCEPT(InvalidParameterError,
                          "radius must be in range [0, PI/2] deg");
    }
    if (radius == 0.0) {
        return 0.0 * radians;
    }
    centerPhi = clampPhi(centerPhi);
    if (std::fabs(centerPhi.asRadians()) + radius.asRadians() > HALFPI - POLE_EPSILON) {
        return PI*(1 + 2.0*DBL_EPSILON) * radians;
    }
    double y = std::sin(radius.asRadians());
    double x = std::sqrt(std::fabs(std::cos(centerPhi.asRadians() - radius.asRadians()) *
                                   std::cos(centerPhi.asRadians() + radius.asRadians())));
    return std::fabs(std::atan(y / x)) * radians;
}

/** Converts from spherical coordinates, proper motions, parallax and
  * radial velocity to position (AU) and velocity (AU/day) 3-vectors.
  *
  * @par
  * See the following paper for details:
  *
  * @par
  * Rigorous Computation of Proper Motions and their Effects on Star Positions
  * Eichhorn, H. & Rust, A.
  * Journal: Astronomische Nachrichten, volume 292, p.37
  * Bibliographic Code: 1970AN....292...37E
  *
  * @par
  * A more accurate approach (including special-relativistic effects) is
  * discussed in Stumpff, P., 1985, Astron.Astrophys. 144, 232-240 and
  * implemented in the IAU SOFA astronomy library iauStarpv C routine
  * available here: http://www.iausofa.org/2009_1231_C/sofa/starpv.c
  *
  * @par
  * Note that the input and output data are for an observer situated at the
  * solar system barycenter.
  *
  * @sa earthPosition(double)
  */
void positionAndVelocity(
    Eigen::Vector3d &p, ///< [out] position, AU
    Eigen::Vector3d &v, ///< [out] velocity, AU/day
    Angle ra,           ///< right ascension
    Angle decl,         ///< declination
    double muRa,        ///< proper motion in RA, rad/day
    double muDecl,      ///< proper motion in Dec, rad/day
    double vRadial,     ///< radial velocity, AU/day
    Angle parallax      ///< parallax
) {
    // distance (AU)
    double r = 1.0 / parallax.asRadians();

    // convert to position and velocity vector (AU, AU/day)
    double sinRa = sin(ra.asRadians());
    double cosRa = cos(ra.asRadians());
    double sinDecl = sin(decl.asRadians());
    double cosDecl = cos(decl.asRadians());
    double s = r*cosDecl;
    double t = r*muDecl*sinDecl;
    double u = cosDecl*vRadial;
    p = Eigen::Vector3d(s*cosRa, s*sinRa, r*sinDecl);
    v = Eigen::Vector3d(cosRa*u - p.y()*muRa - cosRa*t,
                        sinRa*u + p.x()*muRa - sinRa*t,
                        sinDecl*vRadial + s*muDecl);
    if (v.squaredNorm() > 0.25*C_AU_PER_DAY*C_AU_PER_DAY) {
        throw LSST_EXCEPT(RuntimeError,
                          "star velocity vector magnitude exceeds half "
                          "the speed of light");
    }
}

/** Converts the input position vector, which need not have unit magnitude,
  * to spherical coordinates (rad).
  */
Eigen::Vector2d const cartesianToSpherical(Eigen::Vector3d const &v) {
    double d2 = v.x()*v.x() + v.y()*v.y();
    double theta = (d2 == 0.0) ? 0.0 : std::atan2(v.y(), v.x());
    double phi = (v.z() == 0.0) ? 0.0 : std::atan2(v.z(), std::sqrt(d2));
    return Eigen::Vector2d(theta, phi);
}

}}} // namespace lsst::ap::utils

