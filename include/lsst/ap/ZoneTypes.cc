// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of zone related classes.
 *
 * @ingroup ap
 */

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

#include "boost/scoped_array.hpp"

#include "lsst/pex/exceptions.h"

#include "SpatialUtil.h"
#include "ZoneTypes.h"


// -- lsst::ap::ZoneEntry<C> ----------------

template <typename ChunkT>
inline lsst::ap::ZoneEntry<ChunkT>::ZoneEntry(
    double const ra,
    double const dec,
    Data * const data,
    Chunk * const chunk,
    int const index
) :
    _data(data),
    _flags(0),
    _index(index),
    _chunk(chunk)
{
    _ra  = raToScaledInteger(ra);
    _dec = decToScaledInteger(dec);
    double raRad  = radians(ra);
    double decRad = radians(dec);
    double cosDec = std::cos(decRad);
    _x = std::cos(raRad)*cosDec;
    _y = std::sin(raRad)*cosDec;
    _z = std::sin(decRad);
}


// -- lsst::ap::ZoneEntryArray<EntryT> ----------------

template <typename EntryT>
lsst::ap::ZoneEntryArray<EntryT>::ZoneEntryArray()
    : _entries(0), _size(0), _capacity(0), _zone(0), _deltaRa(0) {}


template <typename EntryT>
lsst::ap::ZoneEntryArray<EntryT>::~ZoneEntryArray() {
    if (_entries != 0) {
        std::free(_entries);
        _entries = 0;
    }
}


/** Initializes the zone, allocating space for the given number of entries. */
template <typename EntryT>
void lsst::ap::ZoneEntryArray<EntryT>::init(int const capacity) {
    if (capacity > 0) {
        std::size_t const nb = sizeof(EntryT) * capacity;
        EntryT * entries = static_cast<EntryT *>(std::malloc(nb));
        if (entries == 0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::MemoryException,
                              "failed to allocate zone");
        }
        _entries = entries;
        _capacity = capacity;
    }
}


/** Sorts the zone entries on ra. */
template <typename EntryT>
void lsst::ap::ZoneEntryArray<EntryT>::sort() {
    std::sort(_entries, _entries + _size);
}


/** Increases the size of the underlying array of entries by roughly 25% (and by at least 1). */
template <typename EntryT>
void lsst::ap::ZoneEntryArray<EntryT>::grow() {
    int cap = _capacity >> 2;
    cap = _capacity + (64 > cap ? 64 : cap);
    EntryT * entries  = static_cast<EntryT *>(std::realloc(_entries, sizeof(EntryT)*cap));
    if (entries == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::MemoryException,
                          "failed to increase zone capacity");
    }
    _entries  = entries;
    _capacity = cap;
}


/** Prepares for a distance based match with the given radius */
template <typename EntryT>
void lsst::ap::ZoneEntryArray<EntryT>::computeMatchParams(
    ZoneStripeChunkDecomposition const & zsc,
    double const radius
) {
    double d1 = std::fabs(zsc.getZoneDecMin(_zone));
    double d2 = std::fabs(zsc.getZoneDecMax(_zone));
    _deltaRa  = deltaRaToScaledInteger(maxAlpha(radius, d1 < d2 ? d2 : d1));
}


/**
 * Given a functor that implements @code bool operator()(EntryT const &) @endcode ,
 * removes any entry @a e where @c filter(e) returns @c false from the zone.
 *
 * @return     the number of entries that were removed.
 */
template <typename EntryT>
    template <typename FilterT>
int lsst::ap::ZoneEntryArray<EntryT>::pack(FilterT & filter) {
    int src = 0;
    int dst = 0;

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
    return src - dst;
}


/**
 * Given a functor that implements @code void operator()(EntryT const &) @endcode ,
 * applies it to every entry in the zone.
 */
template <typename EntryT>
    template <typename FunctionT>
void lsst::ap::ZoneEntryArray<EntryT>::apply(FunctionT & function) {
    for (int i = 0; i < _size; ++i) {
        function(_entries[i]);
    }
}


// -- ZoneIndex<EntryT> ----------------

template <typename EntryT>
lsst::ap::ZoneIndex<EntryT>::ZoneIndex(
    int const zonesPerDegree,
    int const zonesPerStripe,
    int const maxEntriesPerZoneEstimate
) :
    lsst::daf::base::Citizen(typeid(*this)),
    _zsc(zonesPerDegree, zonesPerStripe, maxEntriesPerZoneEstimate),
    _zones(),
    _capacity(0),
    _minZone(0),
    _maxZone(-1)
{}


/// Removes all entries from every zone in the index.
template <typename EntryT>
void lsst::ap::ZoneIndex<EntryT>::clear() {
    for (int i = 0; i < _capacity; ++i) {
        _zones[i].clear();
    }
}


/// Returns the number of entries in the index.
template <typename EntryT>
int lsst::ap::ZoneIndex<EntryT>::size() const {
    int sz = 0;
    for (int i = 0; i <= _maxZone - _minZone; ++i) {
        sz += _zones[i].size();
    }
    return sz;
}


/// Sets the range of declination values the index will accept data for.
template <typename EntryT>
void lsst::ap::ZoneIndex<EntryT>::setDecBounds(double const minDec, double const maxDec) {
    int minZone = _zsc.decToZone(minDec);
    int maxZone = _zsc.decToZone(maxDec);
    if (maxZone < minZone) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "min/max zone ids inverted");
    }
    int const cap = maxZone - minZone + 1;
    if (cap >= _capacity) {
        int const worst = _zsc.getMaxEntriesPerZoneEstimate();
        boost::scoped_array<Zone> zones(new Zone[cap]);
        int i = 0;
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

    int i = 0;
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


/// Prepares for distance based matches of the given maximum radius
template <typename EntryT>
void lsst::ap::ZoneIndex<EntryT>::computeMatchParams(double const radius) {
    int const numZones = _maxZone - _minZone + 1;
    for (int t = 0; t < numZones; ++t) {
        _zones[t].computeMatchParams(_zsc, radius);
    }
}


/// Sorts each zone in the index (on right ascension)
template <typename EntryT>
void lsst::ap::ZoneIndex<EntryT>::sort() {
    int const numZones = _maxZone - _minZone + 1;
#if LSST_AP_HAVE_OPEN_MP
#   pragma omp parallel for default(shared) \
               schedule(static,8)
#endif
    for (int t = 0; t < numZones; ++t) {
        _zones[t].sort();
    } // end of parallel for
}


/**
 * Given a functor that implements @code bool operator()(EntryT const &) @endcode ,
 * removes any entry @a e where @c filter(e) returns @c false from the index.
 *
 * @return     the number of entries that were removed.
 */
template <typename EntryT>
    template <typename FilterT>
int lsst::ap::ZoneIndex<EntryT>::pack(FilterT & filter) {
    int const numZones = _maxZone - _minZone + 1;
    int numPacked = 0;
#if LSST_AP_HAVE_OPEN_MP
#   pragma omp parallel for default(shared) \
               reduction(+:numPacked) \
               schedule(static,8)
#endif
    for (int z = 0; z < numZones; ++z) {
        int np = _zones[z].pack<FilterT>(filter);
        numPacked = numPacked + np;
    } // end of parallel for
    return numPacked;
}


/**
 * Calls a functor implementing @code void operator()(EntryT const &) @endcode
 * on every entry in the index.
 */
template <typename EntryT>
    template <typename FunctionT>
void lsst::ap::ZoneIndex<EntryT>::apply(FunctionT & function) {
    int const numZones = _maxZone - _minZone + 1;
    for (int z = 0; z < numZones; ++z) {
        _zones[z].apply<FunctionT>(function);
    }
}


#endif // LSST_AP_ZONE_TYPES_CC
