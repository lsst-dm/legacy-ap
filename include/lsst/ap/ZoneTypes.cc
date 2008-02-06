// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ZoneTypes.cc
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_ZONE_TYPES_CC
#define LSST_AP_ZONE_TYPES_CC

#include <cstdlib>
#include <cmath>

#if LSST_AP_HAVE_OPEN_MP
#   include <omp.h>
#endif

#include <new>
#include <stdexcept>
#include <algorithm>

#include <boost/scoped_array.hpp>

#include "Exceptions.h"
#include "SpatialUtil.h"
#include "ZoneTypes.h"


namespace lsst {
namespace ap {


// -- ZoneEntry<C> ----------------

template <typename C>
ZoneEntry<C>::ZoneEntry(
    DataType  * const data,
    ChunkType * const chunk,
    int32_t     const index
) :
    _data(data),
    _flags(0),
    _index(index),
    _chunk(chunk)
{
    double ra  = data->getRa();
    double dec = data->getDec();
    _ra  = raToScaledInteger(ra);
    _dec = decToScaledInteger(dec);
    ra  = radians(ra);
    dec = radians(dec);
    double cdec = std::cos(dec);
    _x = std::cos(ra)*cdec;
    _y = std::sin(ra)*cdec;
    _z = std::sin(dec);
}


// -- Zone<EntryType> ----------------

template <typename EntryType>
Zone<EntryType>::Zone() : _entries(0), _size(0), _capacity(0), _zone(0), _deltaRa(0) {}


template <typename EntryType>
Zone<EntryType>::~Zone() {
    if (_entries != 0) {
        std::free(_entries);
        _entries = 0;
    }
}


/*! Initializes the zone, allocating space for the given number of entries. */
template <typename EntryType>
void Zone<EntryType>::init(int32_t const capacity) {
    if (capacity > 0) {
        size_t const nb = sizeof(EntryType) * capacity;
        EntryType * entries = static_cast<EntryType *>(std::malloc(nb));
        if (entries == 0) {
            throw std::bad_alloc();
        }
        _entries = entries;
        _capacity = capacity;
    }
}


/*! Sorts the zone entries on ra. */
template <typename EntryType>
void Zone<EntryType>::sort() {
    std::sort(_entries, _entries + _size);
}


/*! Increases the size of the underlying array of entries by roughly 25% (and by at least 1). */
template <typename EntryType>
void Zone<EntryType>::grow() {
    int32_t cap = _capacity >> 2;
    cap = _capacity + (64 > cap ? 64 : cap);
    EntryType * entries  = static_cast<EntryType *>(std::realloc(_entries, sizeof(EntryType)*cap));
    if (entries == 0) {
        throw std::bad_alloc();
    }
    _entries  = entries;
    _capacity = cap;
}


/*! Prepares for a distance based match with the given radius */
template <typename EntryType>
void Zone<EntryType>::computeMatchParams(
    ZoneStripeChunkDecomposition const & zsc,
    double                       const   radius
) {
    double d1 = std::fabs(zsc.getZoneDecMin(_zone));
    double d2 = std::fabs(zsc.getZoneDecMax(_zone));
    _deltaRa  = deltaRaToScaledInteger(maxAlpha(radius, d1 < d2 ? d2 : d1));
}


/*!
    Given a functor that implements \code bool operator()(EntryType const &) \endcode ,
    removes any entry \a e where \c filter(e) returns \c false from the zone.

    \return     the number of entries that were removed.
 */
template <typename EntryType>
    template <typename FilterType>
size_t Zone<EntryType>::pack(FilterType & filter) {
    int32_t src = 0;
    int32_t dst = 0;

    for ( ; src < _size; ++src) {
        if (!filter(_entries[src])) {
            continue;
        }
        if (dst != src) {
            _entries[dst] = _entries[src];
        }
        ++dst;
    }
    _size = dst;
    return static_cast<size_t>(src - dst);
}


/*!
    Given a functor that implements \code void operator()(EntryType const &) \endcode ,
    applies it to every entry in the zone.
 */
template <typename EntryType>
    template <typename FunctionType>
void Zone<EntryType>::apply(FunctionType & function) {
    for (int32_t i = 0; i < _size; ++i) {
        function(_entries[i]);
    }
}


// -- ZoneIndex<EntryType> ----------------

template <typename EntryType>
ZoneIndex<EntryType>::ZoneIndex(
    int32_t const zonesPerDegree,
    int32_t const zonesPerStripe,
    int32_t const maxEntriesPerZoneEstimate
) :
    lsst::mwi::data::Citizen(typeid(*this)),
    _zsc(zonesPerDegree, zonesPerStripe, maxEntriesPerZoneEstimate),
    _zones(),
    _capacity(0),
    _minZone(0),
    _maxZone(-1)
{}


/*! Removes all entries from every zone in the index. */
template <typename EntryType>
void ZoneIndex<EntryType>::clear() {
    for (int32_t i = 0; i < _capacity; ++i) {
        _zones[i].clear();
    }
}


/*! Returns the number of entries in the index. */
template <typename EntryType>
int32_t ZoneIndex<EntryType>::size() const {
    int32_t sz = 0;
    for (int32_t i = 0; i <= _maxZone - _minZone; ++i) {
        sz += _zones[i].size();
    }
    return sz;
}


/*! Sets the range of declination values the index will accept data for. */
template <typename EntryType>
void ZoneIndex<EntryType>::setDecBounds(double const minDec, double const maxDec) {
    int32_t minZone = _zsc.decToZone(minDec);
    int32_t maxZone = _zsc.decToZone(maxDec);
    if (maxZone < minZone) {
        LSST_AP_THROW(InvalidParameter, "min/max zone ids inverted");
    }
    int32_t const cap = maxZone - minZone + 1;
    if (cap >= _capacity) {
        int32_t const worst = _zsc.getMaxEntriesPerZoneEstimate();
        boost::scoped_array<ZoneType> zones(new ZoneType[cap]);
        int32_t i = 0;
        for ( ; i < _capacity; ++i) {
            zones[i] = _zones[i];
        }
        for ( ; i < cap; ++i) {
            zones[i].init(worst);
        }
        // transfer zone ownership from old zone array to new array
        for (i = 0; i < _capacity; ++i) {
            _zones[i]._entries  = 0;
            _zones[i]._size     = 0;
            _zones[i]._capacity = 0;
        }
        using std::swap;
        swap(_zones, zones);
        _capacity = cap;
    }

    int32_t i = 0;
    for ( ; i <= maxZone - minZone; ++i) {
        _zones[i].clear();
        _zones[i]._zone = i + minZone;
    }
    for ( ; i < _capacity; ++i) {
        _zones[i].clear();
    }
    _minZone = minZone;
    _maxZone = maxZone;
}


/*! Prepares for distance based matches of the given maximum radius */
template <typename EntryType>
void ZoneIndex<EntryType>::computeMatchParams(double const radius) {
    int32_t const numZones = _maxZone - _minZone + 1;
    for (int32_t t = 0; t < numZones; ++t) {
        _zones[t].computeMatchParams(_zsc, radius);
    }
}


/*! Sorts each zone in the index (on right ascension) */
template <typename EntryType>
void ZoneIndex<EntryType>::sort() {
    int32_t const numZones = _maxZone - _minZone + 1;
#if LSST_AP_HAVE_OPEN_MP
#   pragma omp parallel for default(shared) \
               schedule(static,8)
#endif
    for (int32_t t = 0; t < numZones; ++t) {
        _zones[t].sort();
    } // end of parallel for
}


/*!
    Given a functor that implements \code bool operator()(EntryType const &) \endcode ,
    removes any entry \a e where \c filter(e) returns \c false from the index.

    \return     the number of entries that were removed.
 */
template <typename EntryType>
    template <typename FilterType>
size_t ZoneIndex<EntryType>::pack(FilterType & filter) {
    size_t        numPacked = 0;
    int32_t const numZones  = _maxZone - _minZone + 1;
#if LSST_AP_HAVE_OPEN_MP
#   pragma omp parallel for default(shared) \
               reduction(+:numPacked) \
               schedule(static,8)
#endif
    for (int32_t z = 0; z < numZones; ++z) {
        size_t np = _zones[z].pack<FilterType>(filter);
        numPacked = numPacked + np;
    } // end of parallel for
    return numPacked;
}


/*!
    Calls a functor implementing \code void operator()(EntryType const &) \endcode
    on every entry in the index.
 */
template <typename EntryType>
    template <typename FunctionType>
void ZoneIndex<EntryType>::apply(FunctionType & function) {
    int32_t const numZones = _maxZone - _minZone + 1;
    for (int32_t z = 0; z < numZones; ++z) {
        _zones[z].apply<FunctionType>(function);
    }
}


}} // end of namespace lsst::ap

#endif // LSST_AP_ZONE_TYPES_CC
