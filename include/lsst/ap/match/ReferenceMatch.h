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

/** @file
  * @brief  Matching against simulated reference catalog positions.
  */
#ifndef LSST_AP_MATCH_REFERENCEMATCH_H
#define LSST_AP_MATCH_REFERENCEMATCH_H

#include <vector>

#include "lsst/afw/geom/Angle.h"

#include "../utils/CsvControl.h"
#include "CatalogControl.h"
#include "ExposureInfo.h"


namespace lsst { namespace ap { namespace match {

// -- Functions to perform the match

/** Matches a declination sorted reference catalog (stored as a CSV file)
  * to a table of positions.
  *
  * Note that a reduction for parallax from barycentric to geocentric place is 
  * applied to reference catalog entries with parallax above parallaxThresh.
  * To disable this reduction, use a large threshold (e.g. +Inf).
  *
  * @param[in] refFile          Declination sorted reference catalog CSV file name.
  * @param[in] refControl       Reference catalog CSV file properties.
  * @param[in] refDialect       CSV dialect of reference catalog CSV file.
  * @param[in] posFile          Declination sorted position CSV file name.
  * @param[in] posControl       Position CSV file properties.
  * @param[in] posDialect       CSV dialect of position CSV file.
  * @param[in] outFile          Output file name.
  * @param[in] outDialect       Output file CSV dialect.
  * @param[in] radius           Match radius.
  * @param[in] parallaxThresh   Parallax threshold.
  * @param[in] outputRefExtras  Output proper-motion/parallax corrected reference
  *                             object position and associated flags in match records?
  * @param[in] truncateOutFile  Truncate outFile before appending to it?
  */
void referenceMatch(
    std::string                 const &refFile,
    CatalogControl              const &refControl,
    lsst::ap::utils::CsvControl const &refDialect,
    std::string                 const &posFile,
    CatalogControl              const &posControl,
    lsst::ap::utils::CsvControl const &posDialect,
    std::string                 const &outFile,
    lsst::ap::utils::CsvControl const &outDialect,
    lsst::afw::geom::Angle      const  radius=2.0*lsst::afw::geom::arcseconds,
    lsst::afw::geom::Angle      const  parallaxThresh=0.01*lsst::afw::geom::arcseconds,
    bool                               outputRefExtras=true,
    bool                               truncateOutFile=false);

/** Computes the number of times a reference catalog should have been observed in
  * each filter with an ideal observatory, given a set of exposures. The per-filter
  * observation counts are appended as columns "<filter>Cov", in order of filter ID.
  * Note that filter IDs are required to be contiguous integers 0, 1, .... N - 1.
  *
  * Reference catalog entries not falling on any of the given exposures are dropped
  * from the output.
  *
  * Note that a reduction for parallax from barycentric to geocentric place is 
  * applied to reference catalog entries with parallax above parallaxThresh.
  * To disable this reduction, use a large threshold (e.g. +Inf).
  *
  * @param[in] exposures        Exposures to filter against - reordered by the call.
  * @param[in] refFile          Declination sorted reference catalog CSV file name.
  * @param[in] refControl       CSV dialect of reference catalog CSV file. 
  * @param[in] refDialect       CSV dialect of reference catalog CSV file.
  * @param[in] outFile          Output file name.
  * @param[in] outDialect       Output file CSV dialect.
  * @param[in] parallaxThresh   Parallax threshold
  * @param[in] truncateOutFile  Truncate outFile before appending to it?
  */
void referenceFilter(
    std::vector<ExposureInfo::Ptr>    &exposures,
    std::string                 const &refFile,
    CatalogControl              const &refControl,
    lsst::ap::utils::CsvControl const &refDialect,
    std::string                 const &outFile,
    lsst::ap::utils::CsvControl const &outDialect,
    lsst::afw::geom::Angle      const  parallaxThresh=0.01*lsst::afw::geom::arcseconds,
    bool                               truncateOutFile=false);

}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_REFERENCEMATCH_H

