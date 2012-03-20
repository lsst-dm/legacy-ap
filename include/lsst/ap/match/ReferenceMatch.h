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
    bool                               truncateOutFile=false);

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

