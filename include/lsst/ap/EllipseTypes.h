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
 

/**
 * @file
 * @brief   Types involved in algorithms for matching points inside ellipses.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_ELLIPSE_TYPES_H
#define LSST_AP_ELLIPSE_TYPES_H

#include <cassert>

#include <vector>

#include "Common.h"
#include "SpatialUtil.h"


namespace lsst { namespace ap {

/**
 * @brief   Contains spatial information for a single ellipse on the unit sphere (sky).
 *
 * A pointer to the actual data object gives access to ancillary fields.
 */
template <typename DataT>
class Ellipse {
public :

    DataT   * _data; ///< pointer to data object (not owned by this object!)
    Ellipse * _next; ///< pointer to next active ellipse in search

    int _minZone;             ///< minimum zone of ellipse bounding-box
    int _maxZone;             ///< maximum zone of ellipse bounding-box
    boost::uint32_t _ra;      ///< right ascension of ellipse center
    boost::uint32_t _deltaRa; ///< width (in right ascension) of ellipse
    boost::int32_t _minDec;   ///< minimum declination of ellipse bounding-box
    boost::int32_t _maxDec;   ///< maximum declination of ellipse bounding-box

    double _sinDec;    ///< sine of ellipse center dec
    double _cosDec;    ///< cosine of ellipse center dec
    double _sinRa;     ///< sine of ellipse center ra
    double _cosRa;     ///< cosine of ellipse center ra
    double _sinPa;     ///< sine of ellipse position angle
    double _cosPa;     ///< cosine of ellipse position angle
    double _invMinor2; ///< 1/(smia*smia), where smia is the ellipse semi-minor axis length (rad)
    double _invMajor2; ///< 1/(smaa*smaa), where smaa is the ellipse semi-major axis length (rad)

    explicit Ellipse(DataT & data);

    /// Returns @c true if the ellipse contains the given unit vector
    bool contains(double const x, double const y, double const z) const {
        // get coords of input point in (N,E) basis
        double xne = _cosDec*z - _sinDec*(_sinRa*y + _cosRa*x);
        double yne = _cosRa*y - _sinRa*x;
        // rotate by negated position angle
        double xr  = _sinPa*yne + _cosPa*xne;
        double yr  = _cosPa*yne - _sinPa*xne;
        // now in position to use standard 2D axis-aligned formulation for an ellipse
        return (xr*xr*_invMajor2 + yr*yr*_invMinor2 <= 1.0);
    }
};

template <typename DataT>
inline void swap(Ellipse<DataT> & a, Ellipse<DataT> & b) {
    a.swap(b);
}

template <typename DataT>
inline bool operator< (Ellipse<DataT> const & a, Ellipse<DataT> const & b) {
    return a._minZone < b._minZone;
}

template <typename DataT>
inline bool operator== (Ellipse<DataT> const & a, Ellipse<DataT> const & b) {
    return a._minZone == b._minZone;
}

template <typename DataT>
inline bool operator< (int32_t const a, Ellipse<DataT> const & b) {
    return a < b._minZone;
}

template <typename DataT>
inline bool operator< (Ellipse<DataT> const & a, int32_t const b) {
    return a._minZone < b;
}

template <typename DataT>
inline bool operator== (int32_t const a, Ellipse<DataT> const & b) {
    return a == b._minZone;
}

template <typename DataT>
inline bool operator== (Ellipse<DataT> const & a, int32_t const b) {
    return a._minZone == b;
}

/** @brief  Comparison functor for Ellipse pointers that orders ellipses by minimum overlapping zone. */
template <typename DataT>
struct EllipsePtrLessThan :
    std::binary_function<Ellipse<DataT> const *, Ellipse<DataT> const *, bool>
{
    bool operator() (Ellipse<DataT> const * a, Ellipse<DataT> const * b) {
        return a->_minZone < b->_minZone;
    }
};


/**
 * @brief   A list of ellipses, implemented using std::vector.
 *
 * Supports the in-ellipse cross matching algorithms.
 */
template <typename DataT>
class EllipseList : public std::vector<Ellipse<DataT> > {
public :
    EllipseList() : std::vector<Ellipse<DataT> >() {}

    /** Creates a list of ellipses from the given data objects. */
    EllipseList(DataT * const begin, DataT * const end);

    void prepareForMatch(ZoneStripeChunkDecomposition const & zsc);
};


}} // end of namespace lsst::ap

#include "EllipseTypes.cc"

#endif // LSST_AP_ELLIPSE_TYPES_H

