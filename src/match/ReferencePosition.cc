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

#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

using lsst::afw::geom::Angle;
using lsst::afw::geom::radians;
using lsst::afw::coord::IcrsCoord;

using lsst::ap::utils::angularSeparation;
using lsst::ap::utils::cartesianToIcrs;
using lsst::ap::utils::clampPhi;
using lsst::ap::utils::maxAlpha;


namespace lsst { namespace ap { namespace match {

/** Clears the motion parameters of this reference position.
  */
void ReferencePosition::clearMotion() {
    _p = _sc.getVector().asEigen();
    _v = Eigen::Vector3d::Zero();
    _parallax = 0.0 * radians;
    _minDecl = _sc.getLatitude();
    _maxDecl = _sc.getLatitude();
    _minRa = _sc.getLongitude();
    _maxRa = _sc.getLongitude();
    _flags = 0;
}

/** Sets the motion parameters of this reference position.
  */
void ReferencePosition::setMotion(
    double muRa,      ///< Rate of change of RA (true or coordinate angle),
                      ///  rad/day
    double muDecl,    ///< Declination rate of change, rad/day 
    Angle parallax,   ///< Parallax, rad
    double vRadial,   ///< Radial velocity, AU/day
    bool trueAngle,   ///< Is muRa dRA/dt*cos(decl) (@c true)
                      ///  or dRA/dt (@c false)?
    bool parallaxCor  ///< Apply parallax corrections in getPosition()?
) {
    double sr = std::sin(_sc.getLongitude().asRadians());
    double cr = std::cos(_sc.getLongitude().asRadians()); 
    double sd = std::sin(_sc.getLatitude().asRadians());
    double cd = std::cos(_sc.getLatitude().asRadians());
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
    double r = 1.0 / parallax.asRadians();
    double s = r*cd;
    double t = r*sd*muDecl;
    double u = cd*vRadial;
    _p.x() = s*cr;
    _p.y() = s*sr;
    _p.z() = r*sd;
    _v.x() = cr*u - _p.y()*muRa - cr*t;
    _v.y() = sr*u + _p.x()*muRa - sr*t;
    _v.z() = sd*vRadial + s*muDecl;
    if (_v.squaredNorm() > 0.25*C_AU_PER_DAY*C_AU_PER_DAY) {
        throw LSST_EXCEPT(pexExcept::RuntimeError,
                          "star velocity vector magnitude exceeds half "
                          "the speed of light");
    }
    _parallax = parallax;
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
        Angle r = std::max(angularSeparation(m, p1),
                           angularSeparation(m, p2));
        if ((_flags & PARALLAX_COR) != 0) {
            r += 2.0*_parallax;
        }
        IcrsCoord sc = cartesianToIcrs(m);
        Angle alpha = maxAlpha(r, sc.getLatitude());
        _minDecl = clampPhi(sc.getLatitude() - r);
        _maxDecl = clampPhi(sc.getLatitude() + r);
        _minRa = sc.getLongitude() - alpha;
        _maxRa = sc.getLongitude() + alpha;
    }
}


double ReferencePosition::getMinCoord0() const {
    return static_cast<double>(_minRa);
}

double ReferencePosition::getMaxCoord0() const {
    return static_cast<double>(_maxRa);
}

double ReferencePosition::getMinCoord1() const {
    return static_cast<double>(_minDecl);
}

double ReferencePosition::getMaxCoord1() const {
    return static_cast<double>(_maxDecl);
}

}}} // namespace lsst::ap::match

