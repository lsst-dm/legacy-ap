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
 * @brief   Useful astronomical constants.
 */

#ifndef LSST_AP_CONSTANTS_H
#define LSST_AP_CONSTANTS_H

namespace lsst { namespace ap { namespace {

double const J2000_MJD     = 51544.5;
double const DAYS_PER_JY   = 365.25;         ///< Days per Julian year
double const METERS_PER_AU = 149597870.7e3;
double const SEC_PER_JD    = 86400.0;        ///< Seconds per Julian day
double const C_AU_PER_DAY  = 173.1446326847; ///< Speed of light in AU/day

}}} // end of namespace lsst::ap::<anonymous>

#endif  // LSST_AP_CONSTANTS_H
