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

// Define alignment macro if the compiler supports it
#if LSST_AP_HAVE_ALIGNED
#   define LSST_AP_ALIGNED(x) __attribute__((aligned(x)))
#else
#   define LSST_AP_ALIGNED(x) 
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

double const J2000_MJD     = 51544.5;
double const DAYS_PER_JY   = 365.25;        ///< Days per Julian year
double const METERS_PER_AU = 149597870.7e3; ///< Meters per AU
double const SEC_PER_JD    = 86400.0;       ///< Seconds per Julian day
double const C_AU_PER_DAY  = 173.1446326847; ///< Speed of light in AU/day

}}} // end of namespace lsst::ap::<anonymous>

#endif  // LSST_AP_COMMON_H
