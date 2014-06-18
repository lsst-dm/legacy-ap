// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2012 LSST Corporation.
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
  * @brief CatalogControl implementation .
  */
#include "lsst/ap/match/CatalogControl.h"

#include <limits>

#include "lsst/pex/exceptions.h"
#include "lsst/utils/ieee.h"


using std::numeric_limits;
using std::string;
using std::vector;
using lsst::pex::exceptions::InvalidParameterError;

namespace lsst { namespace ap { namespace match {

CatalogControl::CatalogControl() :
    minEpoch(numeric_limits<double>::quiet_NaN()),
    maxEpoch(numeric_limits<double>::quiet_NaN()),
    maxParallax(numeric_limits<double>::quiet_NaN()),
    maxAngularVelocity(numeric_limits<double>::quiet_NaN()),
    fieldNames(),
    outputFields(),
    idColumn("id"),
    epochColumn("epoch"),
    epoch(51544.5), // J2000.0 MJD
    raColumn("ra"),
    raScale(0.0174532925199432957692369076849), // rad / deg
    declColumn("decl"),
    declScale(0.0174532925199432957692369076849), // rad / deg
    muRaColumn("muRa"),
    muRaScale(1.32734751843815467101961424328e-11), // (rad / (Julian day)) / (milliarcsec / (Julian yr))
    muRaTrueAngle(true),
    muDeclColumn("muDecl"),
    muDeclScale(1.32734751843815467101961424328e-11), // (rad / (Julian day)) / (milliarcsec / (Julian yr))
    parallaxColumn("parallax"),
    parallaxScale(4.84813681109535993589914102357e-9), // rad / milliarcsec
    vRadialColumn("vRad"),
    vRadialScale(0.000577548327363993690854050371186) // (AU / day) / (km / s)
{
    outputFields.push_back(string("*"));
}

CatalogControl::~CatalogControl() { }

void CatalogControl::validate() const {
    if (lsst::utils::isinf(minEpoch) || lsst::utils::isinf(maxEpoch)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Infinite catalog minEpoch and/or maxEpoch");
    }
    if (lsst::utils::isinf(maxParallax)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Infinite catalog maxParallax");
    }
    if (lsst::utils::isinf(maxAngularVelocity)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Infinite catalog maxAngularVelocity");
    }
    if (!lsst::utils::isfinite(epoch)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Catalog default epoch is not finite");
    }
    if (!lsst::utils::isfinite(raScale) || !lsst::utils::isfinite(declScale)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Catalog raScale and/or declScale (converts from catalog units "
            "to radians) is not finite");
    }
    if (!lsst::utils::isfinite(muRaScale) ||
        !lsst::utils::isfinite(muDeclScale)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Catalog muRaScale and/or muDeclScale (converts from catalog "
            "units to radians per Julian day) is not finite");
    }
    if (!lsst::utils::isfinite(parallaxScale)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Catalog parallaxScale (converts from catalog units to radians) "
            "is not finite");
    }
    if (!lsst::utils::isfinite(vRadialScale)) {
        throw LSST_EXCEPT(InvalidParameterError,
            "Catalog vRadialScale (converts from catalog units to AU per day) "
            "is not finite");
    }
}

}}} // namespace lsst::ap::match

