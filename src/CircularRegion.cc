// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of the CircularRegion class.
 *
 * @ingroup ap
 */

#include "lsst/pex/exceptions.h"

#include "lsst/ap/CircularRegion.h"


namespace ex = lsst::pex::exceptions;

lsst::ap::CircularRegion::CircularRegion(
    double const ra,
    double const dec,
    double const radius
) :
    _centerRa(ra),
    _centerDec(dec),
    _radius(radius)
{
    if (ra < 0.0 || ra >= 360.0) {
        throw LSST_EXCEPT(ex::RangeErrorException,
                          "right ascension must be in range [0, 360) degrees");
    }
    if (dec < -90.0 || dec > 90.0) {
        throw LSST_EXCEPT(ex::RangeErrorException,
                          "declination must be in range  [-90, 90] degrees");
    }
    if (radius < 0.0 || radius > 90.0) {
        throw LSST_EXCEPT(ex::RangeErrorException,
                          "circle radius must be in range  [0, 90] degrees");
    }
}

