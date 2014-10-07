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
  * @brief  Class for simulated reference catalog positions.
  */
#ifndef LSST_AP_MATCH_REFERENCEPOSITION_H
#define LSST_AP_MATCH_REFERENCEPOSITION_H

#include <cstddef>
#include "Eigen/Core"

#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/coord/Coord.h"

#include "../constants.h"
#include "../utils/SpatialUtils.h"
#include "../utils/EarthPosition.h"
#include "BBox.h"


namespace lsst { namespace ap { namespace match {

/** Position related parameters of a simulated reference catalog source.
  * There are no errors - these are inputs to the image simulator.
  */
class ReferencePosition : public BBox {
public:
    enum Flags {
        MOVING       = 0x1, ///< Set if the reference position has proper motion
        PARALLAX     = 0x2, ///< Set if the reference position has
                            ///  parallax \> MIN_PARALLAX
        PARALLAX_COR = 0x4  ///< Set if parallax corrections are applied
                            ///  by getPosition().
    };

    static double const MIN_PARALLAX; // rad

    inline ReferencePosition(int64_t id,
                             lsst::afw::geom::Angle const ra,
                             lsst::afw::geom::Angle const dec,
                             double epoch=J2000_MJD);

    virtual ~ReferencePosition() { }

    void clearMotion();    
    void setMotion(double muRa,
                   double muDecl,
                   lsst::afw::geom::Angle parallax,
                   double vRadial,
                   bool trueAngle,
                   bool parallaxCor);

    void setTimeRange(double epoch1, double epoch2);

    inline int64_t getId() const; 
    inline double getEpoch() const;
    inline int getFlags() const;
    inline lsst::afw::coord::IcrsCoord const & getSphericalCoords() const;
    inline Eigen::Vector3d const & getPosition() const;
    inline Eigen::Vector3d const & getVelocity() const;

    inline Eigen::Vector3d const getPosition(double epoch) const;
    inline Eigen::Vector3d const getPosition(
        double epoch, Eigen::Vector3d const &origin) const;

    // BBox API
    virtual double getMinCoord0() const;
    virtual double getMaxCoord0() const;
    virtual double getMinCoord1() const;
    virtual double getMaxCoord1() const; 

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    lsst::afw::coord::IcrsCoord _sc; ///< (ra, decl) at _epoch, ICRS rad
    int64_t _id;
    double _epoch;       ///< epoch of reference position, MJD
    Eigen::Vector3d _p;  ///< (x, y, z) at _epoch
    Eigen::Vector3d _v;  ///< (dx/dt, dy/dt, dz/dt)
    lsst::afw::geom::Angle _parallax; ///< parallax, rad
    // BBox coordinates
    lsst::afw::geom::Angle _minDecl;
    lsst::afw::geom::Angle _maxDecl;
    lsst::afw::geom::Angle _minRa;
    lsst::afw::geom::Angle _maxRa; 
    int _flags;          ///< Bit-wise OR of Flags
};

}}} // namespace lsst::ap::match

#include "ReferencePosition.cc"

#endif // LSST_AP_MATCH_REFERENCEPOSITION_H

