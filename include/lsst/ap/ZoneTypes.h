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
 * @brief   Classes for zone entries, zones, and zone indexes.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_ZONE_TYPES_H
#define LSST_AP_ZONE_TYPES_H

#include "boost/noncopyable.hpp"
#include "boost/scoped_array.hpp"

#include "lsst/daf/base/Citizen.h"

#include "Common.h"
#include "SpatialUtil.h"


namespace lsst { namespace ap {

/**
 * @brief   Contains spatial information for a single point used during cross-matching.
 *
 * A pointer to the underlying data object gives access to ancillary fields (e.g. colors,
 * magnitudes, etc...).
 */
template <typename ChunkT>
struct LSST_AP_LOCAL ZoneEntry {
    typedef ChunkT Chunk;
    typedef typename ChunkT::Entry Data;

    Data * _data;   ///< Pointer to the corresponding data object
    boost::uint32_t _ra;    ///< scaled right ascension of entity position
    boost::int32_t  _dec;   ///< scaled declination of entity position
    boost::uint32_t _flags; ///< Reserved
    boost::int32_t  _index; ///< Index of the data object in the chunk
    Chunk * _chunk; ///< Pointer to chunk containing the data object
    double  _x;     ///< unit vector x coordinate of entity position
    double  _y;     ///< unit vector y coordinate of entity position
    double  _z;     ///< unit vector z coordinate of entity position

    inline ZoneEntry(
        double const ra,
        double const dec,
        Data * const data,
        Chunk * const chunk,
        int const index
    );
};

template <typename ChunkT>
LSST_AP_LOCAL inline bool operator< (ZoneEntry<ChunkT> const & a, ZoneEntry<ChunkT> const & b) {
    return a._ra < b._ra;
}

template <typename ChunkT>
LSST_AP_LOCAL inline bool operator== (ZoneEntry<ChunkT> const & a, ZoneEntry<ChunkT> const & b) {
    return a._ra == b._ra;
}

template <typename ChunkT>
LSST_AP_LOCAL inline bool operator< (boost::uint32_t const a, ZoneEntry<ChunkT> const & b) {
    return a < b._ra;
}

template <typename ChunkT>
LSST_AP_LOCAL inline bool operator< (ZoneEntry<ChunkT> const & a, boost::uint32_t const b) {
    return a._ra < b;
}

template <typename ChunkT>
LSST_AP_LOCAL inline bool operator== (boost::uint32_t const a, ZoneEntry<ChunkT> const & b) {
    return a == b._ra;
}

template <typename ChunkT>
LSST_AP_LOCAL inline bool operator== (ZoneEntry<ChunkT> const & a, boost::uint32_t const b) {
    return a._ra == b;
}


/**
 * @brief  Stores entries inside a single zone (a narrow declination stripe)
 *         in a sorted array.
 */
template <typename EntryT>
struct LSST_AP_LOCAL ZoneEntryArray {
    typedef typename EntryT::Chunk Chunk;
    typedef typename EntryT::Data  Data;

    EntryT * _entries;
    int _size;
    int _capacity;
    int _zone;
    boost::uint32_t _deltaRa;

    ZoneEntryArray();
    ~ZoneEntryArray();

    void init(int const capacity);

    /** Inserts the given data item into the zone. */
    void insert(double const ra, double const dec, Data * const data, Chunk * const chunk, int const index) {
        int const sz = _size;
        if (sz == _capacity) {
            grow();
        }
        new(&_entries[sz]) EntryT(ra, dec, data, chunk, index);
        _size = sz + 1;
    }

    void sort();

    void grow();

    /** Returns the number of entries in the zone. */
    int size() const { return _size; }

    /** Empties the zone (without deallocating/shrinking memory). */
    void clear() { _size = 0; }

    /** Finds the last entry with ra less than or equal to the specified value. */
    int findLte(boost::uint32_t const ra) {
        EntryT const * const entries = _entries;
        int const last = _size - 1;

        int sz = last + 1;
        int i  = last;

        while (sz > 0) {
            int mid = sz >> 1;
            if (ra < entries[i - mid]) {
                i  -= mid + 1;
                sz -= mid + 1;
            } else {
                sz = mid;
            }
        }
        // if no entry was found, wrap to the largest entry
        return (i < 0) ? last : i;
    }

    /** Finds the first entry with ra greater than or equal to the specified value. */
    int findGte(boost::uint32_t const ra) const {
        EntryT const * const entries = _entries;
        int const end = _size;

        int sz = end;
        int i  = 0;

        while (sz > 0) {
            int mid = sz >> 1;
            if (entries[i + mid] < ra) {
                i  += mid + 1;
                sz -= mid + 1;
            } else {
                sz = mid;
            }
        }
        // if no entry was found, wrap to the smallest entry
        return (i == end) ? 0 : i;
    }

    void computeMatchParams(ZoneStripeChunkDecomposition const & zsc, double const radius);

    template <typename FilterT> int pack (FilterT & filter);
    template <typename FunctionT> void apply(FunctionT & function);
};


/** @brief  Container for a sequence of adjacent zones. */
template <typename EntryT>
class LSST_AP_LOCAL ZoneIndex :
    public  lsst::daf::base::Citizen,
    private boost::noncopyable
{
public :

    typedef typename EntryT::Chunk Chunk;
    typedef typename EntryT::Data  Data;
    typedef ZoneEntryArray<EntryT> Zone;

    ZoneIndex(
        int const zonesPerDegree,
        int const zonesPerStripe,
        int const maxEntriesPerZoneEstimate
    );

    void clear();

    int size() const;

    void setDecBounds(double const minDec, double const maxDec);

    void computeMatchParams(double const radius);

    void sort();

    template <typename FilterT> int pack(FilterT & filter);
    template <typename FunctionT> void apply(FunctionT & function);

    /** Inserts the given data item from the given chunk into the index. */
    void insert(double const ra, double const dec, Data * const data, Chunk * const chunk, int const index) {
        int const zone = _zsc.decToZone(data->getDec());
        if (zone >= _minZone && zone <= _maxZone) {
            _zones[zone - _minZone].insert(ra, dec, data, chunk, index);
        }
    }

    /** Returns the smallest zone id in the index. */
    int getMinZone() const { return _minZone; }

    /** Returns the largest zone id in the index. */
    int getMaxZone() const { return _maxZone; }

    /**
     * Returns a pointer to the zone with the given id,
     * or 0 if the requested zone isn't in the index.
     */
    Zone * getZone(int const zone) {
        if (zone >= _minZone && zone <= _maxZone) {
            return &_zones[zone - _minZone];
        }
        return 0;
    }

    /**
     * Returns a pointer to the first zone in the index within the given id range,
     * or 0 if there is no such zone.
     */
    Zone * firstZone(int const minZone, int const maxZone) {
        if (maxZone < _minZone || minZone > _maxZone) {
            return 0;
        }
        if (minZone <= _minZone) {
            return _zones.get();
        }
        return &_zones[minZone - _minZone];
    }

    /**
     * Returns a pointer to the zone following the last zone in the index within
     * the given id range, or 0 if there is no such zone.
     */
    Zone * endZone(int const minZone, int const maxZone) {
        if (maxZone < _minZone || minZone > _maxZone) {
            return 0;
        }
        if (maxZone >= _maxZone) {
            return &_zones[_maxZone - _minZone + 1];
        }
        return &_zones[maxZone - _minZone + 1];
    }

    ZoneStripeChunkDecomposition const & getDecomposition() const {
        return _zsc;
    }

private :

    ZoneStripeChunkDecomposition _zsc;
    boost::scoped_array<Zone> _zones;
    int _capacity;
    int _minZone;
    int _maxZone;
};


}} // end of namespace lsst::ap

#include "ZoneTypes.cc"

#endif // LSST_AP_ZONE_TYPES_H

