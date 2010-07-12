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
  * @brief Metrics (distance functions) over K-dimensional spaces.
  *
  * A metric is a functor providing operators to computes the distance
  * between two K dimensional points and the minimum distance between
  * two K-dimensional points when given the their k-th coordinate values.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_OPTICS_METRICS_H
#define LSST_AP_CLUSTER_OPTICS_METRICS_H

#include <cmath>

#include "Eigen/Core"

#include "../../Common.h"


namespace lsst { namespace ap { namespace cluster { namespace optics {

/** Metric that computes the squared euclidian distance. Defined
  * over the set of 3-vectors with unit L2 norm.
  */
struct LSST_AP_LOCAL SquaredEuclidianDistanceOverSphere {
    /** Computes the distance between two unit vectors v1, v2.
      */
    double operator()(Eigen::Vector3d const & v1,
                      Eigen::Vector3d const & v2) const
    {
        return (v1 - v2).squaredNorm();
    }

    /** Computes the minimum distance between two unit vectors v1, v2 that
      * have their k-th coordinates fixed to s, t.
      */
    double operator()(double s, double t) const {
        return 2.0 * (1.0 - s * t - std::sqrt((1.0 - s * s) * (1.0 - t * t)));
    }
};

}}}} // lsst::ap::cluster::optics

#endif // LSST_AP_CLUSTER_OPTICS_METRICS_H
