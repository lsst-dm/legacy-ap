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
  * @brief Implementation of the KDTree class template.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_DETAIL_KDTREE_CC
#define LSST_AP_CLUSTER_DETAIL_KDTREE_CC

#include "KDTree.h"

#include <algorithm>
#include <utility>


namespace lsst { namespace ap { namespace cluster { namespace detail {

namespace {

/** @internal
  * Finds the dimension in which the given points have maximum extent.
  * Used to pick a splitting dimension during k-d tree construction.
  *
  * @pre    points != 0
  * @pre    numPoints > 0
  *
  * @return A pair containing the maximum extent of the input points
  *         and the dimension of maximum extent.
  */
template <int K, typename DataT>
std::pair<double, int> maxExtentAndDim(Point<K, DataT> const * const points,
                                       int const numPoints)
{
    typedef Eigen::Matrix<double, K, 1> Vector;

    double minv[K];
    double maxv[K];
    for (int d = 0; d < K; ++d) {
        minv[d] = std::numeric_limits<double>::infinity(); 
        maxv[d] = -std::numeric_limits<double>::infinity();
    }
    for (int i = 0; i < numPoints; ++i) {
        Vector const & v = points[i].coords;
        for (int d = 0; d < K; ++d) {
            minv[d] = std::min(minv[d], v.coeff(d));
            maxv[d] = std::max(maxv[d], v.coeff(d));
        }
    }
    double maxExtent = maxv[0] - minv[0];
    int maxD = 0;
    for (int d = 1; d < K; ++d) {
        double extent = maxv[d] - minv[d];
        if (extent > maxExtent) {
            maxExtent = extent;
            maxD = d;
        }
    }
    return std::make_pair(maxExtent, maxD);
}

/** @internal
  * Orders N-dimensional points along a single dimension.
  */
struct PointCmp {
    int d;

    PointCmp(int dim) : d(dim) { }

    template <int K, typename DataT>
    bool operator()(Point<K, DataT> const & p1, Point<K, DataT> const & p2) const {
        return p1.coords.coeff(d) < p2.coords.coeff(d);
    }
};

} // namespace


/** Creates a new k-d tree over an array of points. The tree construction
  * process modifies the order of points in the array but not the points
  * themselves.
  *
  * @param[in,out] points       Pointer to an array of points.
  * @param[in] numPoints        Number of elements in @c points.
  * @param[in] pointsPerLeaf    Target number of points per leaf node,
  *                             used to determine k-d tree height.
  * @param[in] leafExtentThreshold  If the maximum extent of a k-d tree node
  *                                 along each dimension is below this number,
  *                             then no children are created for the node.
  */
template <int K, typename DataT>
KDTree<K, DataT>::KDTree(Point<K, DataT> * points,
                         int numPoints,
                         int pointsPerLeaf,
                         double leafExtentThreshold) :
    _points(points),
    _numPoints(numPoints),
    _nodes()
{
    if (points == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "pointer to Point array is null");
    }
    if (numPoints <= 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "number of input points must be > 0");
    }
    if (pointsPerLeaf < 1) {
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                          "target number of points per leaf must be > 0");
    }
    // compute tree height
    int h = 0;
    for (; h < MAX_HEIGHT && numPoints / (1 << h) > pointsPerLeaf; ++h) { }
    _height = h;
    int n = (1 << (h + 1)) - 1;
    // allocate tree nodes and build the tree
    boost::scoped_array<KDTreeNode> nodes(new KDTreeNode[n]);
    using std::swap;
    swap(_nodes, nodes);
    build(leafExtentThreshold);
}

template <int K, typename DataT>
KDTree<K, DataT>::~KDTree() { }

/** Locates all points in the k-d tree within distance @c dist of the
  * input query point @c v. Distance computations are performed using
  * @c metric.
  *
  * The result of the query is returned as a single integer index to the
  * first point in range - remaining results are available by traversal of
  * the linked list embedded in the point array. If no points are in range,
  * -1 is returned.
  *
  * @param[in] v        Query point.
  * @param[in] dist     Query distance.
  * @param[in] metric   Distance functor to use for distance computations.
  *
  * @return The index of the first point satsifying the query, or -1 if
  *         none do.
  */
template <int K, typename DataT>
    template <typename MetricT>
int KDTree<K, DataT>::inRange(Vector const & v,
                              double const dist,
                              MetricT const & metric)
{
    bool descend[MAX_HEIGHT];
    int node = 0;
    int h = 0;
    int head = -1;
    int tail = -1;
    while (true) {
        if (_nodes[node].splitDim == -1) {
            // reached a leaf
            int left = 0;
            int right = _nodes[node].right;
            if ((node & (node + 1)) != 0) {
                // node has a left sibling - use it to obtain
                // index of first point in leaf.
                left = _nodes[node - 1].right;
            }
            for (int i = left; i < right; ++i) {
                _points[i].dist = metric(v, _points[i].coords);
            }
            // append results to embedded linked list
            for (int i = left; i < right; ++i) {
                if (_points[i].dist <= dist) {
                    if (tail == -1) {
                        head = i;
                    } else {
                        _points[tail].next = i;
                    }
                    tail = i;
                }
            }
            // move back up the tree
            node = (node - 1) >> 1;
            --h;
            for (; h >= 0 && !descend[h]; --h) {
                node = (node - 1) >> 1;
            }
            if (h < 0) {
                // finished tree traversal
                break;
            }
            descend[h] = false;
            node = (node << 1) + 2;
            ++h;
        } else {
            // determine which children must be visited
            double split = _nodes[node].split;
            double vd = v.coeff(_nodes[node].splitDim);
            if (metric(vd, split) <= dist) {
                // both children must be visited
                descend[h] = true;
                node = (node << 1) + 1;
                ++h;
            } else if (vd < split) {
                // visit left child
                descend[h] = false;
                node = (node << 1) + 1;
                ++h;
            } else {
                // visit right child
                descend[h] = false;
                node = (node << 1) + 2;
                ++h;
            }
        }
    }
    return head;
}

/** @internal
  * Builds k-d tree nodes.
  */
template <int K, typename DataT>
void KDTree<K, DataT>::build(double leafExtentThreshold)
{
    int node = 0;
    int left = 0;
    int right = _numPoints;
    int h = 0;
    while (true) {
        _nodes[node].right = right;
        if (h < _height) {
            // find splitting dimension
            std::pair<double, int> extDim = maxExtentAndDim(_points + left,
                                                            right - left);
            if (extDim.first > leafExtentThreshold) {
                _nodes[node].splitDim = extDim.second;
                // find median of array
                int median = left + ((right - left) >> 1);
                std::nth_element(_points + left, _points + median,
                                 _points + right, PointCmp(extDim.second));
                right = median;
                _nodes[node].split = _points[right].coords.coeff(extDim.second);
                // process left child
                node = (node << 1) + 1;
                ++h;
                continue;
            }
            // node extent is below the subdivision limit: set right index for
            // all right children of node as their left siblings may be valid
            int h2 = h;
            int c = node;
            do {
                c = (c << 1) + 2;
                ++h2;
                _nodes[c].right = right;
            } while (h2 < _height);
        }
        // move up the tree until a left child is found
        left = right;
        for (; h > 0 && (node & 1) == 0; --h) {
            node = (node - 1) >> 1;
        }
        if (h == 0) {
            // tree construction complete!
            break;
        }
        // node is now the index of a left child - process its right sibling
        right = _nodes[(node - 1) >> 1].right;
        node += 1;
    }
}

}}}} // namespace lsst::ap::cluster::detail

#endif // LSST_AP_CLUSTER_DETAIL_KDTREE_CC
