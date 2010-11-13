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
  * @brief  ReferencePosition inline function implementations.
  */
#ifndef LSST_AP_MATCH_REFERENCEPOSITION_CC
#define LSST_AP_MATCH_REFERENCEPOSITION_CC

#include "ReferencePosition.h"


namespace lsst { namespace ap { namespace match {

/** Constructs a stationary, infinitely distant reference position.
  */
ReferencePosition::ReferencePosition(
    int64_t id,   ///< Reference position id
    double ra,    ///< Reference position right ascension, ICRS (radians)
    double decl,  ///< Reference position declination, ICRS (radians)
    double epoch  ///< Reference position epoch, MJD
   ) :
    _sc(ra, decl),
    _id(id),
    _epoch(epoch),
    _p(lsst::ap::utils::sphericalToCartesian(ra, decl)),
    _v(Eigen::Vector3d::Zero()),
    _parallax(0.0),
    _minDecl(decl),
    _maxDecl(decl),
    _minRa(ra),
    _maxRa(ra),
    _flags(0)
{ }

/** Returns the unique integer id of this reference position.
  */
inline int64_t ReferencePosition::getId() const {
    return _id;
}

/** Returns the epoch (MJD) of this reference position.
  */
inline double ReferencePosition::getEpoch() const {
    return _epoch;
}

/** Returns a bit-wise ORed combination of ReferencePosition::Flags describing
  * this reference position.
  */
inline int ReferencePosition::getFlags() const {
   return _flags;
}

/** Returns the ICRS spherical coordinates (rad) of the reference position
  * at epoch getEpoch().
  */
inline Eigen::Vector2d const & ReferencePosition::getSphericalCoords() const {
    return _sc;
}

/** Returns the reference position at epoch getEpoch(). If the
  * reference position flags contain the PARALLAX bit, then this
  * vector is in units of AU and is in the ICRS system (with an origin
  * at the solar-system barycenter). Otherwise, this is a unit-vector -
  * the reference position is treated as infinitely distant from an
  * observer.
  *
  * Coordinates are in the ICRS system for a barycentric observer.
  */
inline Eigen::Vector3d const & ReferencePosition::getPosition() const {
    return _p; 
}

/** Returns the velocity of the reference position. If the
  * reference position flags contain the PARALLAX bit, then this
  * vector is in units of AU/day. Otherwise, units are radians/day
  * (and will not have a radial component).
  */
inline Eigen::Vector3d const & ReferencePosition::getVelocity() const {
    return _v;
}

/** Returns the reference position at epoch @a epoch, accounting for motion
  * and optionally adjusting the coordinates to be geocentric rather than
  * barycentric. The return value is normalized to be a unit 3-vector,
  * which can then be converted to spherical coordinates using
  * lsst::ap::utils::cartesianToSpherical().
  *
  * Note that the change to a geocentric origin is only performed if the
  * the reference position flags contain the PARALLAX bit and @a ssbToGeo
  * is @c true. Otherwise, the reference position is treated as infinitely
  * distant (i.e. the change of origin has no effect),
  */
inline Eigen::Vector3d const ReferencePosition::getPosition(double epoch) const {
    if ((_flags & MOVING) == 0) {
        return _p;
    }
    Eigen::Vector3d p = _p + _v*(epoch - _epoch);
    if ((_flags & PARALLAX_COR) != 0) {
        p -= lsst::ap::utils::earthPosition(epoch);
    }
    return p.normalized();
}

/** Returns the reference position at epoch @a epoch, with a change
  * of origin.
  *
  * @sa getPosition(double) const
  */
inline Eigen::Vector3d const ReferencePosition::getPosition(
    double epoch,
    Eigen::Vector3d const &origin // barycentric coordinate of origin, AU
) const {
    if ((_flags & MOVING) == 0) {
        return _p;
    }
    Eigen::Vector3d p = _p + _v*(epoch - _epoch);
    if ((_flags & PARALLAX_COR) != 0) {
        p -= origin;
    }
    return p.normalized();
}

}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_REFERENCEPOSITION_CC
