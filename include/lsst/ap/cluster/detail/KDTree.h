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
  * @brief Low-level k-d tree class used by the OPTICS implementation.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_DETAIL_KDTREE_H
#define LSST_AP_CLUSTER_DETAIL_KDTREE_H

#include <limits>

#include "boost/scoped_array.hpp"

#include "Eigen/Core"


namespace lsst { namespace ap { namespace cluster { namespace detail {

/** A pointer-less node in a k-d tree. A dimension, splitting value along
  * that dimension, and the index of the point following the last point
  * in the leaf is stored. The index of the first point in the node
  * is obtained from the node to the left at the same level of the tree.
  * By limiting the maximum number of nodes and points to 2^31 - 1,
  * memory usage per node on a typical machine is just 16 bytes - nodes
  * typically will not span cache-lines.
  */
struct KDTreeNode {
    double split;   ///< Splitting value
    int splitDim;   ///< Dimension of splitting value, -1 for leaf nodes
    int right;      ///< Index of first entry to the right of the split

    KDTreeNode() :
        split(std::numeric_limits<double>::quiet_NaN()),
        splitDim(-1),
        right(-1)
    { }
};


/** An entry in the data array to be indexed using a k-d tree.
  * It contains point coordinates, along with the following additional
  * pieces of information:
  *
  *   @li An integer used to embed a singly linked list of range
  *       query results in the data array.
  *   @li A double used to store the distance of the point to the
  *       range query input point.
  *   @li The reachability-distance of the point (defined by the
  *       OPTICS algorithm).
  *   @li An object (or pointer) from which coordinates were obtained.
  *
  * @p
  * Note that if DataT is a pointer, then it is the responsibility
  * of calling code to ensure that the lifetimes of the objects
  * pointed to exceed the lifetimes of any Points referencing them.
  */
template <int K, typename DataT>
struct Point {
    static int const PROCESSED = -2;
    static int const UNPROCESSED = -1;
    
    Eigen::Matrix<double, K, 1> coords; ///< Point coordinates.
    double dist;        ///< Distance to query point.
    double reach;       ///< Reachability distance (for OPTICS).
    DataT data;         ///< Data object.
    int next;           ///< Index of next range query result or -1.
    int state;          ///< State of point ([un]processed or index in seed list)

    Point() :
        coords(),
        dist(std::numeric_limits<double>::quiet_NaN()),
        reach(std::numeric_limits<double>::infinity()),
        data(),
        next(-1),
        state(UNPROCESSED)
    { }

    ~Point() { }
};


/** A pointer-less k-d tree class over an array of k-dimensional Point
  * objects. Points belonging to a node are contiguous in memory.
  * Furthermore, the location of the nodes themselves is implicit:
  * the children of node i are located at positions 2*i + 1 and
  * 2*i + 2 in an underlying array. Nodes therefore need not store
  * pointers to their children, and siblings are contiguous in memory.
  *
  * @p
  * The class supports a simple range query - finding all points within
  * some distance D of a point. The result of this operation is
  * returned as a single integer index to the first Point in range -
  * remaining results are available by traversal of the linked list
  * embedded in the points. Because the results are expected to span a
  * small number of k-d tree leaves and will already have been
  * touched by the range query, the linked list is likely to be
  * cache-resident prior to traversal. However, the consequence of this
  * approach is that a k-d tree and its associated Point array must
  * only be used by a single thread at a time.
  *
  * @p
  * It is also important to note that this class does not own the array
  * of points over which it is defined - it is the responsibility
  * of the caller to ensure that the lifetime of the array exceeds
  * the lifetime of the k-d tree and that the array is not modified
  * while the k-d tree is alive.
  */
template <int K, typename DataT>
class KDTree {
public:
    typedef Eigen::Matrix<double, K, 1> Vector;

    static int const MAX_HEIGHT = 30; ///< Maximum tree height

    KDTree(Point<K, DataT> * points,
           int numPoints,
           int pointsPerLeaf,
           double leafExtentThreshold);
    ~KDTree();

    int size() const {
        return _numPoints;
    }
    int height() const {
        return _height;
    }
    Point<K, DataT> const * getPoints() const {
        return _points;
    }

    template <typename MetricT>
    int inRange(Vector const & v, double const distance, MetricT const & metric);

private:
    Point<K, DataT> * _points;
    int _numPoints;
    int _height;
    boost::scoped_array<KDTreeNode> _nodes;

    void build(double leafExtentThreshold);
};

}}}} // namespace lsst::ap::cluster::detail

#endif // LSST_AP_CLUSTER_DETAIL_KDTREE_H
