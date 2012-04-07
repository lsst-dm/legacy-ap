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
  * @brief  Clustering algorithm control.
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_CLUSTERINGCONTROL_H
#define LSST_AP_CLUSTER_CLUSTERINGCONTROL_H

#include "lsst/pex/config.h"
#include "lsst/afw/geom/Angle.h"


namespace lsst { namespace ap { namespace cluster {

/** @brief Parameters for the clustering algorithm and its internals.
  */
struct ClusteringControl {
    ClusteringControl();
    ~ClusteringControl();

    LSST_CONTROL_FIELD(epsilonArcsec, double,
        "Clustering distance (arcsec) to use when generating clusters with "
        "the OPTICS algorithm.  If a source S has at least minNeighbors "
        "other sources within an angular separation of epsilonArcsec, then "
        "it is always assigned to a cluster.");

    LSST_CONTROL_FIELD(minNeighbors, int,
        "The minimum cardinality of the epsilonArcsec-neighborhood of a source "
        "S for S to be considered a core-source by the OPTICS algorithm. Core "
        "sources are always assigned to a cluster, whereas sources with "
        "epsilonArcsec-neighborhoods containing less than minNeighbors sources "
        "may or may not be.  If such a source is assigned to a cluster, it is "
        "called a border source.  Otherwise, it is called a noise source. "
        "This parameter is essentially a lower bound on the number of times "
        "an astrophysical object must be detected before its detections are "
        "clustered and turned into a catalog entry and can be tuned to avoid "
        "generating too many spurious entries.  To ensure that every source is "
        "assigned to a cluster, set the value to zero.  However, setting the "
        "value to a non-negligeable fraction of the number of times the sky is "
        "covered by the data-set in question will typically result in better "
        "clusters.  Note though that for now, there is no way to adjust the "
        "value to account for data-sets with non uniform coverage.");

    LSST_CONTROL_FIELD(pointsPerLeaf, int,
        "A performance tuning parameter for the k-d tree used internally by "
        "the OPTICS implementation.  The height of the tree is picked such "
        "that no leaf will contain more than pointsPerLeaf sources.  A value "
        "in the tens of sources is generally a good pick.");

    LSST_CONTROL_FIELD(leafExtentThresholdArcsec, double,
        "A performance tuning parameter for the k-d tree used internally by "
        "the OPTICS implementation.  Nodes that have a maximum extent below "
        "this threshold value in each dimension are not subdivided.  The value "
        "should be of the same order as epsilonArcsec - nodes much smaller "
        "than this are useless in the sense that the sources belonging to "
        "such a node become increasingly likely to all lie in the "
        "neighborhood of a query point.  Note also that in the k-d tree "
        "implementation does not store bounding boxes for nodes, meaning "
        "that an entire node cannot be determined to satify a range "
        "query without traversal of its children/contents.  This saves on tree "
        "size, and makes sense for the target use-case because query regions "
        "are typically very small.");

    lsst::afw::geom::Angle const getEpsilon() const {
        return epsilonArcsec * lsst::afw::geom::arcseconds;
    }

    lsst::afw::geom::Angle const getLeafExtentThreshold() const {
        return leafExtentThresholdArcsec * lsst::afw::geom::arcseconds;
    }

    void validate() const;
};

}}} // namespace lsst::ap::cluster

#endif // LSST_AP_CLUSTER_CLUSTERINGCONTROL_H

