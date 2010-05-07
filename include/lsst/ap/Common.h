// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Master header file for the association pipeline.
 *
 * Intended for things of broad utility, e.g. platform normalization and global constants
 *
 * @ingroup ap
 */

#ifndef LSST_AP_COMMON_H
#define LSST_AP_COMMON_H

#include <math.h>
#include <cstddef>

#include "boost/cstdint.hpp"


// -- Shared library support ----------------

#if LSST_AP_HAVE_VISIBILITY
#    define LSST_AP_PUBLIC  __attribute__((visibility("default")))
#    define LSST_AP_LOCAL   __attribute__((visibility("hidden")))
#else
#    define LSST_AP_PUBLIC
#    define LSST_AP_LOCAL
#endif

// Define LSST_AP_API for shared library builds
#if LSST_AP_SHARED_LIBRARY_BUILD
#   define LSST_AP_API LSST_AP_PUBLIC
#else
#   define LSST_AP_API
#endif


namespace lsst { namespace ap { namespace {

/// The radius of an LSST FOV, in degrees.
double const FOV_RADIUS = 1.75;

/**
 * The maximum number of LSST visits in-flight in the association pipeline. In-flight visits
 * are defined as those for which data is actively being read, processed, or written out.
 * @b Must be a power of 2.
 */
int const MAX_VISITS_IN_FLIGHT = 16;

double const DEGREES_PER_RADIAN = 57.2957795130823208767981548141;
double const RADIANS_PER_DEGREE = 0.0174532925199432957692369076849;
double const RADIANS_PER_ARCSEC = 4.84813681109535993589914102357e-6;
double const TWO_PI             = 6.28318530717958647692528676656;
double const ONE_OVER_PI        = 0.318309886183790671537767526745;

}}} // end of namespace lsst::ap::<anonymous>

#endif  // LSST_AP_COMMON_H
