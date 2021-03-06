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
  * @brief  Source processing control.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_SOURCEPROCESSINGCONTROL_H
#define LSST_AP_CLUSTER_SOURCEPROCESSINGCONTROL_H 

#include <string>
#include <vector>
#include "lsst/pex/config.h"


namespace lsst { namespace ap { namespace cluster {

/** @brief Parameters for source processing.
  */
struct SourceProcessingControl {
    SourceProcessingControl();
    ~SourceProcessingControl();

    LSST_CONTROL_FIELD(exposurePrefix, std::string,
        "Prefix for exposure related fields in the output source schema.\n"
        "May be empty.\n");

    LSST_CONTROL_FIELD(clusterPrefix, std::string,
        "Prefix for cluster related fields in the output source schema.\n"
        "May be empty.\n");

    LSST_CONTROL_FIELD(multiBand, bool,
        "If true, sources have been detected/measured on multi-band\n"
        "exposures (e.g. chi-squared coadds) and are not filter\n"
        "specific.\n");

    LSST_CONTROL_FIELD(coadd, bool,
        "If true, sources have been detected/measured on coadd exposures,\n"
        "and have no associated time-stamp/epoch.\n");

    LSST_CONTROL_FIELD(badFlagFields, std::vector<std::string>,
        "A list of flag field names. If an input source has any of the\n"
        "corresponding flag bits set, then the source is considered \"bad\",\n"
        "and does not participate in spatial clustering.\n");

    LSST_CONTROL_FIELD(fluxScale, double,
        "Scaling factor applied to F/F_0 prior to averaging, where F is\n"
        "an uncalibrated source flux and F_0 is the flux of a 0-magnitude\n"
        "object for the corresponding exposure.\n");

    LSST_CONTROL_FIELD(fluxUnit, std::string,
        "Unit of calibrated flux.\n");

    LSST_CONTROL_FIELD(fluxFields, std::vector<std::string>,
        "A list of flux field names which should be carried over from input\n"
        "source tables to output source cluster tables. Input source tables\n"
        "are expected to contain fields '<flux>', '<flux>.err' and\n"
        "'<flux>.flags' for each list entry ('<flux>').\n");

    LSST_CONTROL_FIELD(shapeFields, std::vector<std::string>,
        "A list of shape field names which should be carried over from input\n"
        "source tables to output source cluster tables. Input source tables\n"
        "are expected to contain fields '<shape>', '<shape>.err' and\n"
        "'<shape>.flags' for each list entry ('<shape>').\n");
};

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_SOURCEPROCESSINGCONTROL_H

