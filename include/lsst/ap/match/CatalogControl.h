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
  * @brief  Config parameters for catalogs involved in reference matching.
  */
#ifndef LSST_AP_MATCH_CATALOGCONTROL_H
#define LSST_AP_MATCH_CATALOGCONTROL_H

#include <string>
#include <vector>

#include "lsst/pex/config.h"


namespace lsst { namespace ap { namespace match {


/** Parameters that describe a CSV representation of a catalog involved in
  * reference matching. This can be either the reference catalog, or the
  * catalog of positions being matched against.
  */ 
struct CatalogControl {

    // pre-computed catalog statistics
    LSST_CONTROL_FIELD(minEpoch, double,
        "Minimum epoch of catalog positions, MJD. The default is NaN, "
        "in which case the actual value will be computed (if required) "
        "by scanning the catalog.");

    LSST_CONTROL_FIELD(maxEpoch, double,
        "Maximum epoch of catalog positions, MJD. The default is NaN, "
        "in which case the actual value will be computed (if required) "
        "by scanning the catalog.");              

    LSST_CONTROL_FIELD(maxParallax, double,
        "Maximum parallax (milliarcsec) of any catalog entry. If catalog "
        "entries have no associated parallax, the value should be set to "
        "0.0 or NaN. A value to NaN (the default) will cause the "
        "maximum to be computed (if required) by scanning the catalog.");

    LSST_CONTROL_FIELD(maxAngularVelocity, double,
        "Maximum angular velocity (milliarcsec/yr) of any catalog entry. If "
        "catalog entries have no associated proper motions, the value should "
        "be set to 0.0 or NaN. A value of NaN (the default) will cause the "
        "maximum to be computed (if required) by scanning the catalog. Note "
        "that angular velocity V is related to proper motion components muRa "
        "and muDecl via V = sqrt(muRa**2 + muDecl**2), assuming muRa is in "
        "terms of true angle rather than coordinate angle, i.e. "
        "muRa = dRA/dt*cos(decl).");

    // field names and output options
    LSST_CONTROL_FIELD(fieldNames, std::vector<std::string>,
        "An array of field (column) names for the catalog. If empty, column "
        "names for the catalog will be obtained from the first catalog record.");

    LSST_CONTROL_FIELD(outputFields, std::vector<std::string>,
        "An array of names identifying the catalog fields (columns) to preserve "
        "in output tables. A single value equal to '*' will cause all columns "
        "to be output.");

    // position related column names and scaling factors
    LSST_CONTROL_FIELD(idColumn, std::string,
        "Name of the unique id column for catalog entries. Must be present in "
        "the catalog. NULL values are illegal.");

    LSST_CONTROL_FIELD(epochColumn, std::string,
        "Name of the epoch column for catalog entries. If this name is empty "
        "or not found in the catalog, entries will have their epochs set to "
        "the value of the 'epoch' parameter when read in. NULL epoch values "
        "are illegal.");
    LSST_CONTROL_FIELD(epoch, double,
        "Epoch to assign to catalog entries if the catalog has no epoch "
        "columns. Values must be in MJD, and the default is 51544.5 "
        "(J2000.0 MJD).");

    LSST_CONTROL_FIELD(raColumn, std::string,
        "Name of the ICRS right ascension column. Must be present in the "
        "catalog. NULL values are illegal.");
    LSST_CONTROL_FIELD(raScale, double,
        "Scaling factor to apply to ra column values on reading - must convert "
        "from catalog units to radians. By default, catalog units are assumed "
        "to be degrees.");

    LSST_CONTROL_FIELD(declColumn, std::string,
        "Name of the ICRS declination column. Must be present in the catalog. "
        "NULL values are illegal.");
    LSST_CONTROL_FIELD(declScale, double,
        "Scaling factor to apply to decl column values on reading - must "
        "convert from catalog units to radians. By default, catalog units are "
        "assumed to be degrees.");

    LSST_CONTROL_FIELD(muRaColumn, std::string,
        "Name of the right ascension proper motion column. Need not be present "
        "in the catalog. NULL values are acceptable, and cause the catalog "
        "entry to be treated as stationary and infinitely distant.");
    LSST_CONTROL_FIELD(muRaScale, double,
        "Scaling factor to apply to muRa column values on reading - must "
        "convert from catalog units to radians per Julian day. By default, "
        "catalog units are assumed to be milliarcseconds per Julian year.");
    LSST_CONTROL_FIELD(muRaTrueAngle, bool,
        "True (the default) if the muRa column gives a true angular rate of "
        "change (muRa = dRA/dt*cos(decl)). False if it gives the rate of "
        "change of the coordinate angle (muRa = dRA/dt).");

    LSST_CONTROL_FIELD(muDeclColumn, std::string,
        "Name of the declination proper motion column. Need not be present "
        "in the catalog. NULL values are acceptable, and cause the catalog "
        "entry to be treated as stationary and infinitely distant.");
    LSST_CONTROL_FIELD(muDeclScale, double,
        "Scaling factor to apply to muDecl column values on reading - must "
        "convert from catalog units to radians per Julian day. By default, "
        "catalog units are assumed to be milliarcseconds per Julian year.");

    LSST_CONTROL_FIELD(parallaxColumn, std::string,
        "Name of the parallax column. Need not be present in the catalog. "
        "NULL values are acceptable, and cause the reference catalog entry "
        "to be treated as stationary and infinitely distant.");
    LSST_CONTROL_FIELD(parallaxScale, double,
        "Scaling factor to apply to parallax column values on reading - must "
        "convert from catalog units to radians. By default, catalog units "
        "are assumed to be milliarcseconds.");

    LSST_CONTROL_FIELD(vRadialColumn, std::string,
        "Name of the radial velocity column. Need not be present in the "
        "catalog (in which case zero radial velocity is assumed). NULL "
        "values are also treated as equivalent to 0.0.");
    LSST_CONTROL_FIELD(vRadialScale, double,
        "Scaling factor to apply to vRadial column values on reading - must "
        "convert from catalog units to AU per day. By default, catalog units "
        "are assumed to be km/s.");

    CatalogControl();
    ~CatalogControl();

    void validate() const;
};


}}} // namespace lsst::ap::match

#endif // LSST_AP_MATCH_CATALOGCONTROL_H

