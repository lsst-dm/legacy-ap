// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   ZoneTypes.h
//! \brief  Classes for zone entries, zones, and zone indexes.
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_ZONE_TYPES_H
#define LSST_AP_ZONE_TYPES_H

#include <boost/noncopyable.hpp>
#include <boost/scoped_array.hpp>

#include <lsst/ap/Common.h>
#include <lsst/ap/SpatialUtil.h>


namespace lsst {
namespace ap {


/*!
    \brief  Contains spatial information for a single point used during cross-matching.

    A pointer to the underlying data object gives access to ancillary fields (e.g. colors,
    magnitudes, etc...).
 */
template <typename C>
class ZoneEntry {

public :

    typedef C                     ChunkType;
    typedef typename C::EntryType DataType;

    DataType *  _data;  //!< Pointer to the corresponding data object
    uint32_t    _ra;    //!< scaled right ascension of entity position
    int32_t     _dec;   //!< scaled declination of entity position
    uint32_t    _flags; //!< Reserved
    int32_t     _index; //!< Index of the data object in the chunk
    ChunkType * _chunk; //!< Pointer to chunk containing the data object
    double      _x;     //!< unit vector x coordinate of entity position
    double      _y;     //!< unit vector y coordinate of entity position
    double      _z;     //!< unit vector z coordinate of entity position

    ZoneEntry(DataType * const data, ChunkType * const chunk, int32_t const index);
};

template <typename C>
inline bool operator< (ZoneEntry<C> const & a, ZoneEntry<C> const & b) {
    return a._ra < b._ra;
}

template <typename C>
inline bool operator== (ZoneEntry<C> const & a, ZoneEntry<C> const & b) {
    return a._ra == b._ra;
}

template <typename C>
inline bool operator< (uint32_t const a, ZoneEntry<C> const & b) {
    return a < b._ra;
}

template <typename C>
inline bool operator< (ZoneEntry<C> const & a, uint32_t const b) {
    return a._ra < b;
}

template <typename C>
inline bool operator== (uint32_t const a, ZoneEntry<C> const & b) {
    return a == b._ra;
}

template <typename C>
inline bool operator== (ZoneEntry<C> const & a, uint32_t const b) {
    return a._ra == b;
}


/*! \brief  Contains entries inside a single zone (a narrow declination stripe). */
template <typename EntryType>
class Zone {

public :

    typedef typename EntryType::ChunkType ChunkType;
    typedef typename EntryType::DataType  DataType;

    EntryType * _entries;
    int32_t     _size;
    int32_t     _capacity;
    int32_t     _zone;
    uint32_t    _deltaRa;

    Zone();
    ~Zone();

    void init(int32_t const capacity);

    /*! Inserts the given data item into the zone. */
    void insert(DataType * const data, ChunkType * const chunk, int32_t const index) {
        int32_t const sz = _size;
        if (sz == _capacity) {
            grow();
        }
        new(&_entries[sz]) EntryType(data, chunk, index);
        _size = sz + 1;
    }

    void sort();

    void grow();

    /*! Returns the number of entries in the zone. */
    int32_t size() const { return _size; }

    /*! Empties the zone (without deallocating/shrinking memory). */
    void clear() { _size = 0; }

    /*! Finds the last entry with ra less than or equal to the specified value. */
    int32_t findLte(uint32_t const ra) {
        EntryType const * const entries = _entries;
        int32_t   const         last    = _size - 1;

        int32_t sz = last + 1;
        int32_t i  = last;

        while (sz > 0) {
            int32_t mid = sz >> 1;
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

    /*! Finds the first entry with ra greater than or equal to the specified value. */
    int32_t findGte(uint32_t const ra) const {
        EntryType const * const entries = _entries;
        int32_t   const         end     = _size;

        int32_t sz = end;
        int32_t i  = 0;

        while (sz > 0) {
            int32_t mid = sz >> 1;
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

    template <typename FilterType>   size_t pack (FilterType   & filter);
    template <typename FunctionType> void   apply(FunctionType & function);
};


/*! \brief  Container for a sequence of adjacent zones. */
template <typename EntryType>
class ZoneIndex :
    public  lsst::mwi::data::Citizen,
    private boost::noncopyable
{

public :

    typedef typename EntryType::ChunkType ChunkType;
    typedef typename EntryType::DataType  DataType;
    typedef Zone<EntryType>               ZoneType;

    ZoneIndex(
        int32_t const zonesPerDegree,
        int32_t const zonesPerStripe,
        int32_t const maxEntriesPerZoneEstimate
    );

    void clear();

    int32_t size() const;

    void setDecBounds(double const minDec, double const maxDec);

    void computeMatchParams(double const radius);

    void sort();

    template <typename FilterType>   size_t pack (FilterType   & filter);
    template <typename FunctionType> void   apply(FunctionType & function);

    /*! Inserts the given data item from the given chunk into the index. */
    void insert(DataType * const data, ChunkType * const chunk, int32_t const index) {
        int32_t const zone = _zsc.decToZone(data->getDec());
        if (zone >= _minZone && zone <= _maxZone) {
            _zones[zone - _minZone].insert(data, chunk, index);
        }
    }

    /*! Returns the smallest zone id in the index. */
    int32_t getMinZone() const { return _minZone; }

    /*! Returns the largest zone id in the index. */
    int32_t getMaxZone() const { return _maxZone; }

    /*! Returns a pointer to the zone with the given id, or 0 if the requested zone isn't in the index. */
    ZoneType * getZone(int32_t const zone) {
        if (zone >= _minZone && zone <= _maxZone) {
            return &_zones[zone - _minZone];
        }
        return 0;
    }

    /*! Returns a pointer to the first zone in the index within the given id range,
        or 0 if there is no such zone. */
    ZoneType * firstZone(int32_t const minZone, int32_t const maxZone) {
        if (maxZone < _minZone || minZone > _maxZone) {
            return 0;
        }
        if (minZone <= _minZone) {
            return _zones.get();
        }
        return &_zones[minZone - _minZone];
    }

    /*! Returns a pointer to the zone following the last zone in the index within the given id range,
        or 0 if there is no such zone. */
    ZoneType * endZone(int32_t const minZone, int32_t const maxZone) {
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

    ZoneStripeChunkDecomposition  _zsc;
    boost::scoped_array<ZoneType> _zones;
    int32_t _capacity;
    int32_t _minZone;
    int32_t _maxZone;
};


}} // end of namespace lsst::ap

#include <lsst/ap/ZoneTypes.cc>

#endif // LSST_AP_ZONE_TYPES_H

