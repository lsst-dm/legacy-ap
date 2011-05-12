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
 * @brief   Class and helper functions related to spatial partitioning.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_SPATIAL_UTIL_H
#define LSST_AP_SPATIAL_UTIL_H

#include <cassert>
#include <cmath>

#include <algorithm>
#include <vector>

#include "Common.h"
#include "CircularRegion.h"
#include "RectangularRegion.h"


namespace lsst { namespace ap {

namespace {
    double const RA_DEC_SCALE = 1.19304647111111111111111111111e07; // 1073741824.0/90.0
}

/**
 * @brief   A decomposition of the unit sphere into zones, stripes, and chunks.
 *
 * @em Stripes refer to declination ranges. These are further divided into equal width (in right
 * ascension) @em chunks for the purposes of memory management and/or disk storage. This
 * allows an LSST FOV to be described by a small number (10-100) of spatial chunks that map
 * to some unit of physical storage (e.g. files, database tables).
 *
 * For cross-matching, spatial data is organized into very fine declination ranges, dubbed
 * @em zones, many of which combine to form a stripe. Specifically, a stripe consists of N zones,
 * where N is a positive integer.
 *
 * ZoneStripeChunkDecomposition instances are responsible for mapping declinations/positions
 * to integer zone, stripe, and chunk ids -- these are obtained by simple quantization of
 * declination and right ascension.
 */
class LSST_AP_API ZoneStripeChunkDecomposition {

public :

    ZoneStripeChunkDecomposition(
        int const zonesPerDegree,
        int const zonesPerStripe,
        int const maxEntriesPerZoneEstimate
    );

    ZoneStripeChunkDecomposition(ZoneStripeChunkDecomposition const & zsc);

    ZoneStripeChunkDecomposition & operator=(ZoneStripeChunkDecomposition const & zsc) {
        if (this != &zsc) {
            ZoneStripeChunkDecomposition copy(zsc);
            swap(copy);
        }
        return *this;
    }

    int getNumChunksPerStripe(int const stripeId, double const minWidth) const;

    double stripeAndCircleToRaRange(
        int    const stripeId,
        double const cenRa,
        double const cenDec,
        double const rad
    ) const;

    /**
     * Computes and returns the maximum number of equal-width chunks that can fit into
     * the given declination stripe (where each chunk has a minimum width equal to the
     * stripe height).
     */
    int getNumChunksPerStripe(int const stripeId) const {
        assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of bounds");
        return _chunksPerStripe[stripeId - _minStripe];
    }

    /**
     * Returns the id of the first chunk (the one with the smallest id) in the given
     * declination stripe.
     */
    int getFirstChunkForStripe(int const stripeId) const {
        assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of bounds");
        return stripeId << 16;
    }

    /** Returns the id of the chunk containing the given position */
    int radecToChunk(double const ra, double const dec) const {
        int const s  = decToStripe(dec);
        int const nc = getNumChunksPerStripe(s);
        int c = static_cast<int>(std::floor((ra/360.0)*nc));
        if (c >= nc) {
            --c;
        }
        return (c + getFirstChunkForStripe(s));
    }

    /**
     * Returns the id of the zone that any point with the given declination belongs to.
     * @pre    @code dec >= -90.0 && dec <= 90.0 @endcode
     */
    int decToZone(double const dec) const {
        assert(dec >= -90.0 && dec <= 90.0 && "declination out of bounds");
        int const zoneId = static_cast<int>(std::floor(dec*_zonesPerDegree));
        return (zoneId > _maxZone) ? _maxZone : zoneId;
    }

    /** Returns the minimum declination of points within the given zone. */
    double getZoneDecMin(int const zone) const {
        assert(zone >= _minZone && zone <= _maxZone && "zone id out of bounds");
        double d = static_cast<double>(zone)/_zonesPerDegree;
        return d <= -90.0 ? -90.0 : d;
    }

    /** Returns a non-inclusive upper bound on the declination of points within the given zone. */
    double getZoneDecMax(int const zone) const {
        assert(zone >= _minZone && zone <= _maxZone && "zone id out of bounds");
        double d = static_cast<double>(zone + 1)/_zonesPerDegree;
        return d >= 90.0 ? 90.0 : d;
    }

    /**
     * Returns the id of the stripe that any point with the given declination belongs to.
     * @pre    @code dec >= -90.0 && dec <= 90.0 @endcode
     */
    int decToStripe(double const dec) const {
        assert(dec >= -90.0 && dec <= 90.0 && "declination out of bounds");
        double zoneId = std::floor(dec*_zonesPerDegree);
        return (zoneId > _maxZoneAsDouble) ? _maxStripe :
            static_cast<int>(std::floor(zoneId/_zonesPerStripe));
    }

    /** Returns the id of the stripe the given chunk belongs to. */
    static int chunkToStripe(int const chunkId) {
        return chunkId >> 16;
    }

    /** Returns the sequence number of the given chunk (within its stripe). */
    static int chunkToSequence(int const chunkId) {
        return chunkId & 0x7fff;
    }

    /** Returns the minimum declination of points within the given stripe. */
    double getStripeDecMin(int const stripeId) const {
        assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of bounds");
        double d = static_cast<double>(stripeId*_zonesPerStripe)/_zonesPerDegree;
        return d <= -90.0 ? -90.0 : d;
    }

    /** Returns a non-inclusive upper bound on the declination of points within the given stripe. */
    double getStripeDecMax(int const stripeId) const {
        assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of bounds");
        double d = static_cast<double>((stripeId + 1)*_zonesPerStripe)/_zonesPerDegree;
        return d >= 90.0 ? 90.0 : d;
    }

    /** Returns the smallest zone within the given stripe. */
    int getStripeZoneMin(int const stripeId) const {
        assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of bounds");
        return stripeId*_zonesPerStripe;
    }

    /** Returns the largest zone within the given stripe. */
    int getStripeZoneMax(int const stripeId) const {
        assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of bounds");
        return stripeId*_zonesPerStripe + _zonesPerStripe - 1;
    }

    int getZonesPerStripe() const { return _zonesPerStripe; }

    /**
     * Returns an estimate of the largest number of entries that can fall within the
     * intersection of a zone and an LSST FOV. Always returns a non-negative integer.
     */
    int getMaxEntriesPerZoneEstimate() const {
        return _maxEntriesPerZoneEstimate;
    }

    void swap(ZoneStripeChunkDecomposition & zsc);

private :

    std::vector<int> _chunksPerStripe;

    double _zonesPerDegree;
    double _maxZoneAsDouble;
    int _zonesPerStripe;
    int _maxEntriesPerZoneEstimate;
    int _minZone;
    int _maxZone;
    int _minStripe;
    int _maxStripe;
};

inline void swap(ZoneStripeChunkDecomposition & a, ZoneStripeChunkDecomposition & b) {
    a.swap(b);
}


LSST_AP_API double alpha(double const theta, double const centerDec, double const dec);

LSST_AP_API double maxAlpha(double const theta, double const centerDec);

LSST_AP_API void computeChunkIds(
    std::vector<int>                   & chunkIds,
    CircularRegion               const & region,
    ZoneStripeChunkDecomposition const & zsc,
    int                          const   workerId    = 0,
    int                          const   numWorkers  = 1
);

LSST_AP_API void computeChunkIds(
    std::vector<int>                   & chunkIds,
    RectangularRegion            const & region,
    ZoneStripeChunkDecomposition const & zsc,
    int                          const   workerId    = 0,
    int                          const   numWorkers  = 1
);

/** Converts a right ascension (in degrees) to an integer in the range [0, 2^32) */
inline boost::uint32_t raToScaledInteger(double const ra) {
    assert(ra >= 0.0 && ra < 360.0);
    double d = std::floor(ra*RA_DEC_SCALE);
    return d >= 4294967295.0 ? 4294967295u : static_cast<uint32_t>(d);
}

/** Converts a right ascension delta (in degrees) to an integer */
inline boost::uint32_t deltaRaToScaledInteger(double const delta) {
    assert(delta >= 0.0 && delta <= 180.0);
    return static_cast<int32_t>(std::ceil(delta*RA_DEC_SCALE));
}

/** Converts a declination (in degrees) to an integer between -2^30 and 2^30 */
inline boost::int32_t decToScaledInteger(double const dec) {
    assert(dec >= -90.0 && dec <= 90.0);
    return static_cast<int32_t>(std::floor(dec*RA_DEC_SCALE));
}

/** Converts a declination delta (in degrees) to an integer */
inline boost::int32_t deltaDecToScaledInteger(double const delta) {
    assert(delta >= 0.0 && delta <= 90.0);
    return static_cast<int32_t>(std::ceil(delta*RA_DEC_SCALE));
}

/** Clamps the given declination value to [-90,90]. */
inline double clampDec(double const dec) {
    if (dec <= -90.0) {
        return -90.0;
    } else if (dec >= 90.0) {
        return 90.0;
    }
    return dec;
}


}} // end of namespace lsst::ap

#endif // LSST_AP_SPATIAL_UTIL_H
