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
  * @brief  ClusteringControl implementation.
  * @author Serge Monkewitz
  */
#include "lsst/ap/cluster/ClusteringControl.h"

#include "lsst/pex/exceptions.h"


using lsst::pex::exceptions::InvalidParameterError;

namespace lsst { namespace ap { namespace cluster {

ClusteringControl::ClusteringControl() :
    epsilonArcsec(0.75),
    minNeighbors(2),
    pointsPerLeaf(32),
    leafExtentThresholdArcsec(2.0)
{
    validate();
}

ClusteringControl::~ClusteringControl() { }

void ClusteringControl::validate() const {
    if (epsilonArcsec < 0.0 || epsilonArcsec > 36000.0) {
        throw LSST_EXCEPT(InvalidParameterError,
                          "epsilonArcsec must lie in the range [0, 36000]");
    }
    if (minNeighbors < 0) {
        throw LSST_EXCEPT(InvalidParameterError,
                          "minNeighbors must be non-negative");

    }
    if (pointsPerLeaf <= 0) {
        throw LSST_EXCEPT(InvalidParameterError,
                          "pointsPerLeaf must be positive");
    }
}

}}} // namespace lsst::ap::cluster
