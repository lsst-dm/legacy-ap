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

#include "lsst/pex/policy/Policy.h"
#include "ExposureInfo.h"


namespace lsst { namespace ap { namespace match {

// -- Functions to perform the match

LSST_AP_API void referenceMatch(
    std::string const &refInPath,
    std::string const &posInPath,
    std::string const &matchOutPath,
    lsst::pex::policy::Policy::Ptr refInPolicy=lsst::pex::policy::Policy::Ptr(),
    lsst::pex::policy::Policy::Ptr posInPolicy=lsst::pex::policy::Policy::Ptr(),
    lsst::pex::policy::Policy::Ptr matchPolicy=lsst::pex::policy::Policy::Ptr(),
    bool truncate=false);

LSST_AP_API void referenceFilter(
    std::string const &refInPath,
    std::string const &filtOutPath,
    std::vector<ExposureInfo::Ptr> &exposures,
    lsst::pex::policy::Policy::Ptr refInPolicy=lsst::pex::policy::Policy::Ptr(),
    lsst::pex::policy::Policy::Ptr outPolicy=lsst::pex::policy::Policy::Ptr(),
    bool truncate=false);

}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_REFERENCEMATCH_H

