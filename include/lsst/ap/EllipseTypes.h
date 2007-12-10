// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   EllipseTypes.h
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_ELLIPSE_TYPES_H
#define LSST_AP_ELLIPSE_TYPES_H

#include <cassert>

#include <vector>

#include <lsst/ap/Common.h>
#include <lsst/ap/SpatialUtil.h>


namespace lsst {
namespace ap {


/*!
    \brief  Contains spatial information for a single ellipse on the unit sphere (sky).

    A pointer to the actual data object gives access to ancillary fields.
 */
template <typename DataType>
class Ellipse {

public :

    DataType * _data;        //!< pointer to data object
    Ellipse *  _next;        //!< pointer to next active ellipse in search

    int32_t    _minZone;     //!< minimum zone of ellipse bounding-box
    int32_t    _maxZone;     //!< maximum zone of ellipse bounding-box
    uint32_t   _ra;          //!< right ascension of ellipse center
    uint32_t   _deltaRa;     //!< width (in right ascension) of ellipse
    int32_t    _minDec;      //!< minimum declination of ellipse bounding-box
    int32_t    _maxDec;      //!< maximum declination of ellipse bounding-box

    double     _sinDec;      //!< sine of ellipse center dec
    double     _cosDec;      //!< cosine of ellipse center dec
    double     _sinRa;       //!< sine of ellipse center ra
    double     _cosRa;       //!< cosine of ellipse center ra
    double     _sinPa;       //!< sine of ellipse position angle
    double     _cosPa;       //!< cosine of ellipse position angle
    double     _invMinor2;   //!< 1/(smia*smia), where smia is the ellipse semi-minor axis length (rad)
    double     _invMajor2;   //!< 1/(smaa*smaa), where smaa is the ellipse semi-major axis length (rad)

    explicit Ellipse(DataType & data);

    /*! Returns \c true if the ellipse contains the given unit vector */
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

template <typename DataType>
inline void swap(Ellipse<DataType> & a, Ellipse<DataType> & b) {
    a.swap(b);
}

template <typename DataType>
inline bool operator< (Ellipse<DataType> const & a, Ellipse<DataType> const & b) {
    return a._minZone < b._minZone;
}

template <typename DataType>
inline bool operator== (Ellipse<DataType> const & a, Ellipse<DataType> const & b) {
    return a._minZone == b._minZone;
}

template <typename DataType>
inline bool operator< (int32_t const a, Ellipse<DataType> const & b) {
    return a < b._minZone;
}

template <typename DataType>
inline bool operator< (Ellipse<DataType> const & a, int32_t const b) {
    return a._minZone < b;
}

template <typename DataType>
inline bool operator== (int32_t const a, Ellipse<DataType> const & b) {
    return a == b._minZone;
}

template <typename DataType>
inline bool operator== (Ellipse<DataType> const & a, int32_t const b) {
    return a._minZone == b;
}

/*! \brief  Comparison functor for Ellipse pointers that orders ellipses by minimum overlapping zone. */
template <typename DataType>
struct EllipsePtrLessThan :
    std::binary_function<Ellipse<DataType> const *, Ellipse<DataType> const *, bool>
{
    bool operator() (Ellipse<DataType> const * a, Ellipse<DataType> const * b) {
        return a->_minZone < b->_minZone;
    }
};


/*!
    \brief  A list of ellipses, implemented using std::vector.

    Supports the in-ellipse cross matching algorithms.
 */
template <typename DataType>
class EllipseList {

public :

    typedef Ellipse<DataType> EllipseType;

    typedef typename std::vector<EllipseType>::iterator       iterator;
    typedef typename std::vector<EllipseType>::const_iterator const_iterator;
    typedef typename std::vector<EllipseType>::size_type      size_type;

    EllipseList() : _ellipses() {}

    /*! Creates a list of ellipses from the given data objects. */
    EllipseList(
        DataType * const begin,
        DataType * const end
    );

    EllipseList(EllipseList const & list);

    EllipseList & operator=(EllipseList const & list);

    size_type size()     const { return _ellipses.size();     }
    size_type capacity() const { return _ellipses.capacity(); }
    bool      empty()    const { return _ellipses.empty();    }

    EllipseType const & operator[](size_type const i) const { return _ellipses[i]; }
    EllipseType const * begin() const { return &_ellipses.front();    }
    EllipseType const * end()   const { return &_ellipses.back() + 1; }

    EllipseType & operator[](size_type const i) { return _ellipses[i]; }
    EllipseType * begin() { return &_ellipses.front();    }
    EllipseType * end()   { return &_ellipses.back() + 1; }

    void push_back(DataType & data) { _ellipses.push_back(EllipseType(data)); }
    void pop_back()                 { _ellipses.pop_back();                   }
    void clear()                    { _ellipses.clear();                      }
    void reserve(size_type const n) { _ellipses.reserve(n);                   }
    void swap(EllipseList & list)   { std::swap(_ellipses, list._ellipses);   }

    void prepareForMatch(ZoneStripeChunkDecomposition const & zsc);

private :

    std::vector<EllipseType> _ellipses;
};


}} // end of namespace lsst::ap

#include <lsst/ap/EllipseTypes.cc>

#endif // LSST_AP_ELLIPSE_TYPES_H

