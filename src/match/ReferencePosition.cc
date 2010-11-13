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
  * @brief  Implementation for ReferencePosition class.
  */
#include "lsst/ap/match/ReferencePosition.h"

#include <algorithm>
#include <cmath>

using lsst::ap::utils::angularSeparation;
using lsst::ap::utils::cartesianToSpherical;
using lsst::ap::utils::maxAlpha;
using lsst::ap::utils::sphericalToCartesian;


namespace lsst { namespace ap { namespace match {

/** Clears the motion parameters of this reference position.
  */
void ReferencePosition::clearMotion() {
    _p = sphericalToCartesian(_sc.x(), _sc.y());
    _v = Eigen::Vector3d::Zero();
    _parallax = 0.0;
    _minDecl = _sc.y();
    _maxDecl = _sc.y();
    _minRa = _sc.x();
    _maxRa = _sc.x();
    _flags = 0;
}

/** Sets the motion parameters of this reference position.
  */
void ReferencePosition::setMotion(
    double muRa,      ///< Rate of change of RA (true or coordinate angle),
                      ///  rad/day
    double muDecl,    ///< Declination rate of change, rad/day 
    double parallax,  ///< Parallax, rad
    double vRadial,   ///< Radial velocity, AU/day
    bool trueAngle,   ///< Is muRa dRA/dt*cos(decl) (@c true)
                      ///  or dRA/dt (@c false)?
    bool parallaxCor  ///< Apply parallax corrections in getPosition()?
) {
    double sr = std::sin(_sc.x());
    double cr = std::cos(_sc.x()); 
    double sd = std::sin(_sc.y());
    double cd = std::cos(_sc.y());
    if (trueAngle) {
        muRa = (cd == 0.0) ? 0.0 : muRa/cd;
    }
    // if parallax is tiny, treat this as a position on the celestial sphere
    if (parallax < MIN_PARALLAX) {
        _p.x() = cd*cr;
        _p.y() = cd*sr;
        _p.z() = sd;
        double t = sd*muDecl;
        _v.x() = -_p.y()*muRa - cr*t;
        _v.y() =  _p.x()*muRa - sr*t;
        _v.z() = cd*muDecl;
        _flags = MOVING;
        return;
    }
    double r = 1.0/parallax;
    double s = r*cd;
    double t = r*sd*muDecl;
    _p.x() = s*cr;
    _p.y() = s*sr;
    _p.z() = r*sd;
    _v.x() = _p.x()*vRadial - _p.y()*muRa - cr*t;
    _v.y() = _p.y()*vRadial + _p.x()*muRa - sr*t;
    _v.z() = _p.z()*vRadial + s*muDecl;
    _flags = MOVING | PARALLAX | (parallaxCor ? PARALLAX_COR : 0);
}

/** Sets the bounding box (in spherical coordinates) of the reference position
  * to the bounding box of its path over the given time range. If SSB to
  * geocentric corrections are enabled, the box is additionally padded by twice
  * the parallax.
  *
  * The input epochs need not be ordered.
  */
void ReferencePosition::setTimeRange(double epoch1, double epoch2) {
    if ((_flags & MOVING) != 0) {
        Eigen::Vector3d p1 = _p + _v*(epoch1 - _epoch);
        Eigen::Vector3d p2 = _p + _v*(epoch2 - _epoch);
        Eigen::Vector3d m = p1 + p2;
        double r = std::max(angularSeparation(m, p1),
                            angularSeparation(m, p2));
        if ((_flags & PARALLAX_COR) != 0) {
            r += 2.0*_parallax;
        }
        Eigen::Vector2d sc = cartesianToSpherical(m);
        double alpha = maxAlpha(r, sc.y());
        _minDecl = sc.y() - r;
        _maxDecl = sc.y() + r;
        _minRa = sc.x() - alpha;
        _maxRa = sc.x() + alpha;
    }
}


double ReferencePosition::getMinCoord0() const {
    return _minRa;
}

double ReferencePosition::getMaxCoord0() const {
    return _maxRa;
}

double ReferencePosition::getMinCoord1() const {
    return _minDecl;
}

double ReferencePosition::getMaxCoord1() const {
    return _maxDecl;
}

}}} // namespace lsst::ap::match

