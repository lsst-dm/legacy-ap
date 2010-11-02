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
#include "lsst/ap/util/SpatialUtils.h"

#include <cmath>

#include "lsst/pex/exceptions.h"


namespace lsst { namespace ap { namespace util {

/** Reduces a theta (longitude/right-ascension) range. The resulting
  * range will have <tt> min \> max </tt> if it wraps across the 0/360 degree
  * discontinuity. Valid inputs are:
  *
  * @li  any <tt> min, max </tt> with <tt> min \<= max </tt>
  * @li  <tt> min \> max </tt> so long as
  *      <tt> min \<= 360.0 && max \>= 0.0 </tt>
  */
void thetaRangeReduce(double &min, double &max) {
    if (min > max) {
        if (max < 0.0 || min >= 360.0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                              "Invalid longitude interval");
        }
    } else if (max - min >= 360.0) {
        min = 0.0;
        max = 360.0;
    } else if (min < 0.0 || max >= 360.0) {
        // range reduce
        min = std::fmod(min, 360.0);
        max = std::fmod(max, 360.0);
        if (min < 0.0) { min += 360.0; }
        if (max < 0.0) { max += 360.0; }
    }
}

/** Computes the extent in longitude angle <tt> [-alpha, alpha] </tt> of the
  * circle with radius @a radius and center at latitude angle @a centerPhi.
  *
  * @pre  <tt> radius > 0.0 && radius <= 90.0 </tt>
  * @pre  <tt> centerPhi >= -90.0 && centerPhi <= 90.0 </tt>
  */
double maxAlpha(double radius,   ///< circle radius (deg)
                double centerPhi ///< latitude angle of circle center (deg)
               )
{
    static const double POLE_EPSILON = 1e-6;

    if (radius < 0.0 || radius > 90.0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "radius must be in range [0, 90] deg");
    }
    if (radius == 0.0) {
        return 0.0;
    }
    centerPhi = clampPhi(centerPhi);
    if (std::fabs(centerPhi) + radius > 90.0 - POLE_EPSILON) {
        return 180.0;
    }
    double r = radians(radius);
    double c = radians(centerPhi);
    double y = std::sin(r);
    double x = std::sqrt(std::fabs(std::cos(c - r) * std::cos(c + r)));
    return degrees(std::fabs(std::atan(y / x)));
}

/** Applies the classical space-motion transform to obtain the barycentric
  * coordinates of the input position at @a toEpoch.
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
  * solar system barycenter - before matching these positions against the
  * source positions computed by the LSST pipelines, a reduction for parallax
  * from barycentric to topocentric place should be applied.
  *
  * @return  The barycentric coordinates of the star (AU) at the desired epoch.
  *
  * @sa earthPosition(double)
  */
lsst::afw::geom::Point3D const starProperMotion(
    double ra,          ///< right ascension at @a fromEpoch (ICRS), rad
    double decl,        ///< declination at @a fromEpoch (ICRS), rad
    double muRa,        ///< proper motion in RA, rad per Julian day
    double muDecl,      ///< proper motion in Dec, rad per Julian day
    double vRadial,     ///< radial velocity, AU per Julian day
    double parallax,    ///< parallax, rad
    double fromEpoch,   ///< starting epoch, JD or MJD
    double toEpoch      ///< ending epoch, JD or MJD
) {
    // distance (AU)
    double r = 1.0 / parallax;

    // convert to position and velocity vector (AU, AU/day)
    double sinRa = sin(ra);
    double cosRa = cos(ra);
    double sinDecl = sin(decl);
    double cosDecl = cos(decl);
    double s = r*cosDecl;
    double t = r*muDecl*sinDecl;
    double dt = toEpoch - fromEpoch;
    Eigen::Vector3d p(s*cosRa,
                      s*sinRa,
                      r*sinDecl);
    Eigen::Vector3d v(p.x()*vRadial - p.y()*muRa - cosRa*t,
                      p.y()*vRadial + p.x()*muRa - sinRa*t,
                      p.z()*vRadial              + s*muDecl);
    // compute position at toEpoch
    return lsst::afw::geom::Point3D(p + v*dt);
}

/** Converts the input position vector, which need not have unit magnitude,
  * to spherical coordinates (rad).
  */
lsst::afw::geom::Point2D const cartesianToSpherical(lsst::afw::geom::Point3D const &v) {
    double x = v.coeffRef(0);
    double y = v.coeffRef(1);
    double z = v.coeffRef(2);
    double d2 = x*x + y*y;
    double theta = (d2 == 0.0) ? 0.0 : std::atan2(y, x);
    double phi = (z == 0.0) ? 0.0 : std::atan2(z, std::sqrt(d2));
    return lsst::afw::geom::makePointD(theta, phi);
}

}}} // namespace lsst::ap::util

