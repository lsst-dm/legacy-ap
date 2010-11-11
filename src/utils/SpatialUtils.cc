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

namespace pexExcept = lsst::pex::exceptions;


namespace lsst { namespace ap { namespace utils {

/** Reduces a theta (longitude/right-ascension) range. The resulting
  * range will have <tt> min \> max </tt> if it wraps across the 0/2*M_PI
  * radiandiscontinuity. Valid inputs are:
  *
  * @li  any <tt> min, max </tt> with <tt> min \<= max </tt>
  * @li  <tt> min \> max </tt> so long as
  *      <tt> min \<= 2*M_PI && max \>= 0.0 </tt>
  */
void thetaRangeReduce(double &min, double &max) {
    if (min > max) {
        if (max < 0.0 || min >= 2.0*M_PI) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                              "Invalid longitude angle interval");
        }
    } else if (max - min >= 2.0*M_PI) {
        min = 0.0;
        max = 2.0*M_PI;
    } else if (min < 0.0 || max >= 2.0*M_PI) {
        // range reduce
        min = std::fmod(min, 2.0*M_PI);
        max = std::fmod(max, 2.0*M_PI);
        if (min < 0.0) { min += 2.0*M_PI; }
        if (max < 0.0) { max += 2.0*M_PI; }
    }
}

/** Computes the extent in longitude angle <tt> [-alpha, alpha] </tt> (rad) of
  * the circle with radius @a radius and center at latitude angle @a centerPhi.
  *
  * @pre  <tt> radius >= 0.0 && radius <= M_PI/2 </tt>
  *
  * Note that @a centerPhi is clamped to lie in [-M_PI/2, M_PI/2].
  */
double maxAlpha(double radius,   ///< circle radius (rad)
                double centerPhi ///< latitude angle of circle center (rad)
               )
{
    static const double POLE_EPSILON = 1e-7;

    if (radius < 0.0 || radius > M_PI*0.5) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException,
                          "radius must be in range [0, M_PI/2] deg");
    }
    if (radius == 0.0) {
        return 0.0;
    }
    centerPhi = clampPhi(centerPhi);
    if (std::fabs(centerPhi) + radius > M_PI*0.5 - POLE_EPSILON) {
        return M_PI*(1 + 2.0*DBL_EPSILON);
    }
    double y = std::sin(radius);
    double x = std::sqrt(std::fabs(std::cos(centerPhi - radius) *
                                   std::cos(centerPhi + radius)));
    return std::fabs(std::atan(y / x));
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
void positionAndVelocity(Eigen::Vector3d &p, ///< [out] position, AU
                         Eigen::Vector3d &v, ///< [out] velocity, AU/day
                         double ra,          ///< right ascension, rad
                         double decl,        ///< declination, rad
                         double muRa,        ///< proper motion in RA, rad/day
                         double muDecl,      ///< proper motion in Dec, rad/day
                         double vRadial,     ///< radial velocity, AU/day
                         double parallax     ///< parallax, rad
                        )
{
    // distance (AU)
    double r = 1.0 / parallax;

    // convert to position and velocity vector (AU, AU/day)
    double sinRa = sin(ra);
    double cosRa = cos(ra);
    double sinDecl = sin(decl);
    double cosDecl = cos(decl);
    double s = r*cosDecl;
    double t = r*muDecl*sinDecl;
    p = Eigen::Vector3d(s*cosRa, s*sinRa, r*sinDecl);
    v = Eigen::Vector3d(p.x()*vRadial - p.y()*muRa - cosRa*t,
                        p.y()*vRadial + p.x()*muRa - sinRa*t,
                        p.z()*vRadial              + s*muDecl);
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
