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
 * @brief   Implementation of Ellipse and EllipseList classes.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_ELLIPSE_TYPES_CC
#define LSST_AP_ELLIPSE_TYPES_CC

#include <cmath>

#include "boost/scoped_array.hpp"

#include "EllipseTypes.h"
#include "SpatialUtil.h"


// -- Ellipse<DataT> ----------------

template <typename DataT>
lsst::ap::Ellipse<DataT>::Ellipse(DataT & data) {
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
lsst::ap::EllipseList<DataT>::EllipseList(DataT * const begin, DataT * const end)
    : std::vector<Ellipse<DataT> >()
{
    std::ptrdiff_t const n = end - begin;
    for (std::ptrdiff_t i = 0; i < n; ++i) {
        push_back(Ellipse<DataT>(begin[i]));
    }
}

/**
 * Prepares the list of ellipses for matching against a list of positions --
 * finds zone bounds for each ellipse and sorts them in order of minimum zone.
 */
template <typename DataT>
void lsst::ap::EllipseList<DataT>::prepareForMatch(ZoneStripeChunkDecomposition const & zsc) {
    typedef typename std::vector<DataT>::size_type size_type; 
    size_type const sz = this->size();

    // find min and max zone for each ellipse
    for (size_type i = 0; i < sz; ++i) {
        Ellipse<DataT> & e = this->operator[](i);
        double dec  = e._data->getDec();
        double smaa = e._data->getSemiMajorAxisLength();
        // Note - this actually computes the min/max zone of the ellipses bounding circle
        double d = dec - smaa;
        e._minZone = zsc.decToZone(d <= -90.0 ? -90.0 : d);
        d = dec + smaa;
        e._maxZone = zsc.decToZone(d >= 90.0 ? 90.0 : d);
        e._next = 0;
    }

    if (sz > 1) {
        // indirect sort on min zone
        boost::scoped_array<Ellipse<DataT> *> pointers(new Ellipse<DataT> *[sz]);
        for (size_type i = 0; i < sz; ++i) {
            pointers[i] = &this->operator[](i);
        }
        std::sort(pointers.get(), pointers.get() + sz, EllipsePtrLessThan<DataT>());
        // copy ellipses into new array in sorted order
        std::vector<Ellipse<DataT> > ellipses;
        ellipses.reserve(sz);
        for (size_type i = 0; i < sz; ++i) {
            ellipses.push_back(*pointers[i]);
        }
        // replace old array with new sorted one
        this->swap(ellipses);
    }
}

#endif // LSST_AP_ELLIPSE_TYPES_CC
