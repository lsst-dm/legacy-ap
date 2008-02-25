// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of Ellipse and EllipseList classes.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_ELLIPSE_TYPES_CC
#define LSST_AP_ELLIPSE_TYPES_CC

#include <cmath>

#include <stdexcept>

#include <boost/scoped_array.hpp>

#include "EllipseTypes.h"
#include "SpatialUtil.h"


namespace lsst {
namespace ap {


// -- Ellipse<DataT> ----------------

template <typename DataT>
Ellipse<DataT>::Ellipse(DataT & data) {
    double ra   = data.getRa();
    double dec  = data.getDec();
    double pa   = data.getPositionAngle();
    double smia = data.getSemiMinorAxisLength()/3600.0;
    double smaa = data.getSemiMajorAxisLength()/3600.0;
    _data       = &data;
    _next       = 0;
    _ra         = raToScaledInteger(ra);
    _deltaRa    = deltaRaToScaledInteger(maxAlpha(smaa, dec));
    // Note - this actually computes the min/max dec of the ellipses bounding circle
    double d    = dec - smaa;
    _minDec     = decToScaledInteger(d <= -90.0 ? -90.0 : d);
    d           = dec + smaa;
    _maxDec     = decToScaledInteger(d >= 90.0 ? 90.0 : d);
    ra          = radians(ra);
    dec         = radians(dec);
    pa          = radians(pa);
    _sinDec     = std::sin(dec);
    _cosDec     = std::cos(dec);
    _sinRa      = std::sin(ra);
    _cosRa      = std::cos(ra);
    _sinPa      = std::sin(pa);
    _cosPa      = std::cos(pa);
    smia        = radians(smia);
    smaa        = radians(smaa);
    _invMinor2  = 1.0/(smia*smia);
    _invMajor2  = 1.0/(smaa*smaa);
}


// -- EllipseList<DataT> ----------------

template <typename DataT>
EllipseList<DataT>::EllipseList(DataT * const begin, DataT * const end) : _ellipses() {
    ptrdiff_t const n = end - begin;
    for (ptrdiff_t i = 0; i < n; ++i) {
        push_back(begin[i]);
    }
}


template <typename DataT>
EllipseList<DataT>::EllipseList(EllipseList const & list) : _ellipses(list._ellipses) {}


template <typename DataT>
EllipseList<DataT> & EllipseList<DataT>::operator=(EllipseList const & list) {
    if (this != &list) {
        EllipseList copy(list);
        swap(copy);
    }
    return *this;
}


/**
 * Prepares the list of ellipses for matching against a list of positions --
 * finds zone bounds for each ellipse and sorts them in order of minimum zone.
 */
template <typename DataT>
void EllipseList<DataT>::prepareForMatch(ZoneStripeChunkDecomposition const & zsc) {

    size_type const sz = size();

    // find min and max zone for each ellipse
    for (size_type i = 0; i < sz; ++i) {
        Ellipse & e    = _ellipses[i];
        double    dec  = e._data->getDec();
        double    smaa = e._data->getSemiMajorAxisLength();
        // Note - this actually computes the min/max zone of the ellipses bounding circle
        double d = dec - smaa;
        e._minZone = zsc.decToZone(d <= -90.0 ? -90.0 : d);
        d = dec + smaa;
        e._maxZone = zsc.decToZone(d >= 90.0 ? 90.0 : d);
        e._next    = 0;
    }

    if (sz > 1) {
        // indirect sort on min zone
        boost::scoped_array<Ellipse *> pointers(new Ellipse *[sz]);
        for (size_type i = 0; i < sz; ++i) {
            pointers[i] = &_ellipses[i];
        }
        std::sort(pointers.get(), pointers.get() + sz, EllipsePtrLessThan<DataT>());
        // copy ellipses into new array in sorted order
        std::vector<Ellipse> ellipses;
        ellipses.reserve(sz);
        for (size_type i = 0; i < sz; ++i) {
            ellipses.push_back(*pointers[i]);
        }
        // replace old array with new sorted one
        using std::swap;
        swap(_ellipses, ellipses);
    }
}


}} // end of namespace lsst::ap

#endif // LSST_AP_ELLIPSE_TYPES_CC
