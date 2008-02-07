// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Implementation of the CircularRegion class.
 *
 * @ingroup associate
 */

#include <stdexcept>

#include <lsst/ap/CircularRegion.h>
#include <lsst/ap/Exceptions.h>


namespace lsst {
namespace ap {


CircularRegion::CircularRegion(
    double const ra,
    double const dec,
    double const radius
) :
    _centerRa(ra),
    _centerDec(dec),
    _radius(radius)
{
    if (ra < 0.0 || ra >= 360.0) {
        LSST_AP_THROW(OutOfRange, "right ascension must be in range [0, 360) degrees");
    }
    if (dec < -90.0 || dec > 90.0) {
        LSST_AP_THROW(OutOfRange, "declination must be in range  [-90, 90] degrees");
    }
    if (radius < 0.0 || radius > 90.0) {
        LSST_AP_THROW(OutOfRange, "circle radius must be in range  [0, 90] degrees");
    }
}


}}  // end of namespace lsst::ap

