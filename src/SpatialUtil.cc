// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of spatial utility classes and functions.
 *
 * @ingroup associate
 */

#include <cassert>
#include <stdexcept>
#include <algorithm>

#include <lsst/ap/Exceptions.h>
#include <lsst/ap/SpatialUtil.h>


namespace lsst {
namespace ap {


// -- ZoneStripeChunkDecomposition ----------------

ZoneStripeChunkDecomposition::ZoneStripeChunkDecomposition(
    int32_t const zonesPerDegree,
    int32_t const zonesPerStripe,
    int32_t const maxEntriesPerZoneEstimate
) :
    _chunksPerStripe(),
    _zonesPerDegree(zonesPerDegree),
    _zonesPerStripe(zonesPerStripe),
    _maxEntriesPerZoneEstimate(maxEntriesPerZoneEstimate)
{
    if (zonesPerDegree < 1 || zonesPerDegree > 3600) {
        LSST_AP_THROW(OutOfRange, "zone height must be between 1 arc-second and 1 degree");
    }
    if (zonesPerStripe < 1) {
        LSST_AP_THROW(OutOfRange, "there must be at least 1 zone per stripe");
    }
    if (maxEntriesPerZoneEstimate < 1) {
        LSST_AP_THROW(OutOfRange, "the max entries per zone estimate must be at least 1");
    }
    _minZone         = -90*zonesPerDegree;
    _maxZone         = 90*zonesPerDegree - 1;
    _maxZoneAsDouble = static_cast<double>(_maxZone);

    // C89/C++ standard says: rounding direction of division when either operand is negative
    // is implementation defined. Therefore, implement round to -Infinity by hand.
    int32_t const quo = (-_minZone) / zonesPerStripe;
    int32_t const rem = (-_minZone) % zonesPerStripe;
    _minStripe        = (rem == 0) ? -quo : -1 - quo;
    _maxStripe        = _maxZone/zonesPerStripe;

    int32_t const numStripes = _maxStripe - _minStripe + 1;
    double  const minWidth   = static_cast<double>(zonesPerStripe)/_zonesPerDegree;
    if (numStripes >= 32768) {
        LSST_AP_THROW(OutOfRange, "Requested spatial parameters result in more than 32767 stripes");
    }
    if (getNumChunksPerStripe(0, minWidth) >= 32768) {
        LSST_AP_THROW(OutOfRange, "Requested spatial parameters result in more than 32767 chunks per stripe");
    }

    _chunksPerStripe.reserve(numStripes);
    for (int32_t i = 0; i < numStripes; ++i) {
        _chunksPerStripe.push_back(getNumChunksPerStripe(i + _minStripe, minWidth));
    }
}


ZoneStripeChunkDecomposition::ZoneStripeChunkDecomposition(ZoneStripeChunkDecomposition const & zsc) :
    _chunksPerStripe(zsc._chunksPerStripe),
    _zonesPerDegree (zsc._zonesPerDegree),
    _maxZoneAsDouble(zsc._maxZoneAsDouble),
    _zonesPerStripe (zsc._zonesPerStripe),
    _maxEntriesPerZoneEstimate(zsc._maxEntriesPerZoneEstimate),
    _minZone        (zsc._minZone),
    _maxZone        (zsc._maxZone),
    _minStripe      (zsc._minStripe),
    _maxStripe      (zsc._maxStripe)
{}


void ZoneStripeChunkDecomposition::swap(ZoneStripeChunkDecomposition & zsc) {
    using std::swap;

    if (this != &zsc) {
        swap(_chunksPerStripe, zsc._chunksPerStripe);
        swap(_zonesPerDegree,  zsc._zonesPerDegree);
        swap(_maxZoneAsDouble, zsc._maxZoneAsDouble);
        swap(_maxEntriesPerZoneEstimate, zsc._maxEntriesPerZoneEstimate);
        swap(_zonesPerStripe,  zsc._zonesPerStripe);
        swap(_minZone,         zsc._minZone);
        swap(_maxZone,         zsc._maxZone);
        swap(_minStripe,       zsc._minStripe);
        swap(_maxStripe,       zsc._maxStripe);
   }
}


/**
 * Computes and returns the maximum number of equal-width chunks that can fit
 * into the given declination stripe, where each chunk has the given minimum width.
 * The minimum width for a chunk is defined as the minimum allowable distance between
 * two points in non-adjacent chunks belonging to the same stripe.
 *
 * @param[in] stripeId  the stripe to determine a chunk count for
 * @param[in] minWidth  the minimum width of a chunk (in degrees)
 */
int32_t ZoneStripeChunkDecomposition::getNumChunksPerStripe(
    int32_t const stripeId,
    double const  minWidth
) const {
    assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of range");
    assert(minWidth > 0.0 && minWidth < 90.0 && "chunk width out of range");

    double d1 = std::fabs(getStripeDecMin(stripeId));
    double d2 = std::fabs(getStripeDecMax(stripeId));
    double maxAbsDec = d1 < d2 ? d2 : d1;
    if (maxAbsDec > 89.9) {
        return 1;
    }
    double cosWidth = std::cos(radians(minWidth));
    double sinDec   = std::sin(radians(maxAbsDec));
    double cosDec   = std::cos(radians(maxAbsDec));
    cosWidth = (cosWidth - sinDec*sinDec)/(cosDec*cosDec);
    if (cosWidth < 0) {
        return 1;
    }
    return static_cast<int32_t>(std::floor(TWO_PI/std::acos(cosWidth)));
}


/**
 * Computes the intersection of the given declination stripe and circular region,
 * then returns the range of right ascension values in the intersection as a
 * difference in ra from the circle center.
 * 
 * @param[in]  stripeId id of the declination stripe
 * @param[in]  cenRa    right ascension of circle center (degrees)
 * @param[in]  cenDec   declination of circle center (degrees)
 * @param[in]  rad      radius of circle (degrees)
 *
 * @return  alpha, such that points in the intersection of the input circle and stripe
 *          have right ascensions between [cenRa - alpha, cenRa + alpha].
 */
double ZoneStripeChunkDecomposition::stripeAndCircleToRaRange(
    int32_t const stripeId,
    double  const cenRa,
    double  const cenDec,
    double  const rad
) const {
    assert(stripeId >= _minStripe && stripeId <= _maxStripe && "stripe id out of range");

    // find dec extents of stripe and circle
    double const stripeDecMin = getStripeDecMin(stripeId);
    double const stripeDecMax = getStripeDecMax(stripeId);
    double circleDecMin = cenDec - rad;
    double circleDecMax = cenDec + rad;
    double a            = 0.0;

    if (circleDecMax <= stripeDecMin || circleDecMin >= stripeDecMax) {
        // the circle doesn't intersect the stripe
        a = 0.0;
    } else if (circleDecMax >= 89.9) {
        // the circle contains (or is very close to) the north pole
        if (circleDecMax > 90.0) {
            circleDecMax = 180.0 - circleDecMax;
        }
        if (stripeDecMax >= circleDecMax) {
            a = 180.0;
        } else {
            a = alpha(rad, cenDec, stripeDecMax);
        }
    } else if (circleDecMin <= -89.9) {
        // the circle contains (or is very close to) the south pole
        if (circleDecMin < -90.0) {
            circleDecMin = -180.0 - circleDecMin;
        }
        if (stripeDecMin <= circleDecMin) {
            a = 180.0;
        } else {
            a = alpha(rad, cenDec, stripeDecMin);
        }
    } else {
        // the circle does not contain a pole
        double decOfMaxAlpha = degrees(std::asin(std::sin(radians(cenDec))/std::cos(radians(rad))));
        if (decOfMaxAlpha < stripeDecMin) {
            a = alpha(rad, cenDec, stripeDecMin);
        } else if (decOfMaxAlpha > stripeDecMax) {
            a = alpha(rad, cenDec, stripeDecMax);
        } else {
            // dec of max alpha value is within the stripe
            a = maxAlpha(rad, cenDec);
        }
    }
    return a;
}


// -- Helper functions ----------------

/**
 * Intersects the plane given by z = sin(@a dec) with the circle of radius @a theta and center
 * (0, @a centerDec) on the unit sphere, then returns the right ascension alpha giving the two
 * resulting points: (-alpha, @a dec) and (alpha, @a dec).
 *
 * @pre    @code theta > 0.0 && theta < 10.0 @endcode
 * @pre    @code centerDec >= -90.0 && centerDec <= 90.0 @endcode
 * @pre    @code dec >= -90.0 && dec <= 90.0 @endcode
 *
 * @param[in] theta     the radius of the input circle
 * @param[in] centerDec the declination of the circle center
 * @param[in] dec       the declination of the horizontal plane to
 *                      intersect with the input circle
 *
 * @return  alpha, the largest right ascension of the two points in the intersection
 *          of the input circle with the input plane.
 */
LSST_AP_API double alpha(
    double const theta,
    double const centerDec,
    double const dec
) {
    assert(theta > 0.0 && theta < 10.0 && "radius out of range");
    assert(centerDec >= -90.0 && centerDec <= 90.0 && "declination out of range");
    assert(dec >= -90.0 && dec <= 90.0 && "declination out of range");

    if (std::fabs(centerDec) + theta > 89.9) {
        return 180.0;
    }
    double x = std::cos(radians(theta)) - std::sin(radians(centerDec))*std::sin(radians(dec));
    double u = std::cos(radians(centerDec))*std::cos(radians(dec));
    double y = std::sqrt(std::fabs(u*u - x*x));
    return degrees(std::fabs(std::atan2(y,x)));
}


/**
 * Computes the extent in right ascension [-alpha, alpha] of the circle
 * with radius @a theta and center (0, @a centerDec) on the unit sphere.
 *
 * @pre    @code theta > 0.0 && theta < 10.0 @endcode
 * @pre    @code centerDec >= -90.0 && centerDec <= 90.0 @endcode
 *
 * @param[in] theta     the radius of the circle to find ra extents for
 * @param[in] centerDec the declination of the circle center (in degrees)
 *
 * @return  the largest right ascension of any point on the input circle
 */
LSST_AP_API double maxAlpha(
    double const theta,
    double const centerDec
) {
    assert(theta > 0.0 && theta < 10.0 && "radius out of range");
    assert(centerDec >= -90.0 && centerDec <= 90.0 && "declination out of range");

    if (std::fabs(centerDec) + theta > 89.9) {
        return 180.0;
    }
    double y = std::sin(radians(theta));
    double x = std::sqrt(std::fabs(cos(radians(centerDec - theta))*cos(radians(centerDec + theta))));
    return degrees(std::fabs(std::atan(y/x)));
}


/**
 * Computes identifiers for all chunks in the given ZoneStripeChunkDecomposition that overlap
 * the given region and belong to the specified worker. Chunks belonging to a stripe @c s such
 * that @code s % numWorkers == workerId @endcode belong to the worker identified by @a workerId.
 *
 * @param[out] chunkIds     The list in which to store the computed chunk identifiers.
 * @param[in]  region       The region for which overlapping chunks are to be computed.
 * @param[in]  zsc          A decomposition of the unit sphere into stripes, chunks, and zones.
 * @param[in]  workerId     The integer id of the current worker (in a set of @a numWorkers parallel
 *                          workers).
 * @param[in]  numWorkers   The number of parallel workers.
 */
LSST_AP_API void computeChunkIds(
    std::vector<int64_t>               & chunkIds,
    CircularRegion               const & region,
    ZoneStripeChunkDecomposition const & zsc,
    int                          const   workerId,
    int                          const   numWorkers
) {
    if (numWorkers < 1) {
        LSST_AP_THROW(InvalidParameter, "number of workers must be positive");
    }
    if (workerId < 0 || workerId >= numWorkers) {
        LSST_AP_THROW(
            InvalidParameter,
            "Worker id must be between 0 and N - 1, where N is the total number of workers"
        );
    }

    for (int32_t s = zsc.decToStripe(region.getMinDec()); s <= zsc.decToStripe(region.getMaxDec()); ++s) {

        // round-robin stripes to workers
        int32_t rem = s % numWorkers;
        if (rem < 0) {
            rem += numWorkers;
        }
        if (rem != workerId) {
            continue;
        }

        int32_t const fc = zsc.getFirstChunkForStripe(s);
        int32_t const nc = zsc.getNumChunksPerStripe(s);
        double a = zsc.stripeAndCircleToRaRange(s,
            region.getCenterRa(), region.getCenterDec(), region.getRadius());
        double raMin = region.getCenterRa() - a;
        double raMax = region.getCenterRa() + a;
        bool   wrap  = false;

        if (raMin < 0.0) {
            raMin += 360.0;
            wrap   = true;
        }
        if (raMax >= 360.0) {
            raMax -= 360.0;
            wrap   = true;
        }
        int32_t maxChunk = static_cast<int32_t>(std::floor((raMax/360.0)*nc));
        if (maxChunk == nc) {
            --maxChunk;
        }
        int32_t chunk = static_cast<int32_t>(std::floor((raMin/360.0)*nc));
        if (chunk == nc) {
            --chunk;
        }
        chunk    += fc;
        maxChunk += fc;

        if (raMax < raMin || (raMax == raMin && wrap)) {
            if (chunk == maxChunk) {
                --maxChunk; // avoid adding the same chunk twice
            }
            for ( ; chunk < fc + nc; ++chunk) {
                chunkIds.push_back(chunk);
            }
            for (chunk = 0; chunk <= maxChunk; ++chunk) {
                chunkIds.push_back(chunk);
            }
        } else {
            for ( ; chunk <= maxChunk; ++chunk) {
                chunkIds.push_back(chunk);
            }
        }
    }
}


/**
 * Computes identifiers for all chunks in the given ZoneStripeChunkDecomposition that overlap
 * the given region and belong to the specified worker. Chunks belonging to a stripe @c s such
 * that @code s % numWorkers == workerId @endcode belong to the worker identified by @a workerId.
 *
 * @param[out] chunkIds     The list in which to store the computed chunk identifiers.
 * @param[in]  region       The region for which overlapping chunks are to be computed.
 * @param[in]  zsc          A decomposition of the unit sphere into stripes, chunks, and zones.
 * @param[in]  workerId     The integer id of the current worker (in a set of @a numWorkers parallel
 *                          workers).
 * @param[in]  numWorkers   The number of parallel workers.
 */
LSST_AP_API void computeChunkIds(
    std::vector<int64_t>               & chunkIds,
    RectangularRegion            const & region,
    ZoneStripeChunkDecomposition const & zsc,
    int32_t                      const   workerId,
    int32_t                      const   numWorkers
) {
    if (numWorkers < 1) {
        LSST_AP_THROW(InvalidParameter, "number of workers must be positive");
    }
    if (workerId < 0 || workerId >= numWorkers) {
        LSST_AP_THROW(InvalidParameter,
                      "Worker id must be between 0 and N - 1, where N is the total number of workers");
    }

    double const raMin = region.getMinRa();
    double const raMax = region.getMaxRa();

    for (int32_t s = zsc.decToStripe(region.getMinDec()); s <= zsc.decToStripe(region.getMaxDec()); ++s) {

        // round-robin stripes to workers
        int32_t rem = s % numWorkers;
        if (rem < 0) {
            rem += numWorkers;
        }
        if (rem != workerId) {
            continue;
        }

        int32_t const fc = zsc.getFirstChunkForStripe(s);
        int32_t const nc = zsc.getNumChunksPerStripe(s);

        int32_t maxChunk = fc + static_cast<int32_t>(std::floor((raMax*nc)/360.0)) % nc;
        int32_t chunk    = fc + static_cast<int32_t>(std::floor((raMin*nc)/360.0)) % nc;

        if (raMax < raMin) {
            if (chunk == maxChunk) {
                --maxChunk; // avoid adding the same chunk twice
            }
            for ( ; chunk < fc + nc; ++chunk) {
                chunkIds.push_back(chunk);
            }
            for (chunk = 0; chunk <= maxChunk; ++chunk) {
                chunkIds.push_back(chunk);
            }
        } else {
            for (; chunk <= maxChunk; ++chunk) {
                chunkIds.push_back(chunk);
            }
        }
    }
}


}} // end of namespace lsst::ap
