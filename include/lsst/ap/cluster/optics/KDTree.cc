// -*- lsst-c++ -*-
/** @file
  * @brief Implementation of the KDTree class template.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_CLUSTER_OPTICS_KDTREE_CC
#define LSST_AP_CLUSTER_OPTICS_KDTREE_CC

#include "KDTree.h"

#include <utility>


namespace lsst { namespace ap { namespace cluster { namespace optics {

namespace except = ::lsst::pex::exceptions;

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
  * Partitions the given array of points around the i-th point in the array.
  *
  * @pre    points != 0
  * @pre    numPoints > 0
  * @pre    dim >= 0 && dim < K
  * @pre    i >= 0 && i < numPoints
  *
  * @return The new index of the pivot point.
  */
template <int K, typename DataT>
int partition(Point<K, DataT> * const points,
              int const numPoints,
              int const dim,
              int const i)
{
    Point<K, DataT> pivot = points[i];
    double pivotValue = pivot.coords.coeff(dim);
    points[i] = points[numPoints - 1];
    int j = 0;
    for (int k = 0; k < numPoints - 1; ++k) {
        if (points[k].coords.coeff(dim) < pivotValue) {
            std::swap(points[j], points[k]);
            ++j;
        }
    }
    points[numPoints - 1] = points[j];
    points[j] = pivot;
    return j;
}

///@internal Finds the median of 2 points along dimension @c dim.
template <int K, typename DataT>
inline int median2(Point<K, DataT> * const points, int const dim) {
    return (points[0].coords.coeff(dim) < points[1].coords.coeff(dim)) ? 0 : 1;
}

///@internal Finds the median of 3 points along dimension @c dim.
template <int K, typename DataT>
inline int median3(Point<K, DataT> * const points, int const dim) {
    double v0 = points[0].coords.coeff(dim);
    double v1 = points[1].coords.coeff(dim);
    double v2 = points[2].coords.coeff(dim);
    if (v0 < v1) {
        return v1 < v2 ? 1 : (v0 < v2 ? 2 : 0);
    } else {
        return v1 < v2 ? (v0 < v2 ? 0 : 2) : 1;
    }
}

///@name Median finding lookup tables
//@{
/** @internal
  * Lookup tables for the 4 and 5 element median finding algorithms.
  * Computed with the following python 2.6 script:
  *
  * @code
  * import itertools
  * 
  * def computeLut(n):
  *     nbits = (n * (n - 1)) / 2
  *     lut = [-1] * 2**nbits
  *     array = range(n)
  *     median = array[len(array) >> 1]
  *     for p in itertools.permutations(array):
  *         res = []
  *         for i in xrange(n - 1):
  *             for j in xrange(i + 1, n):
  *                 res.append(1 if p[i] < p[j] else 0)
  *         index = 0
  *         for i in xrange(len(res)):
  *             index += res[i] << (len(res) - i - 1)
  *         lut[index] = p.index(median)
  *     return lut
  * @endcode
  */
signed char const lut4[64] = {
     1, 1,-1, 3, 2,-1, 2, 3,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1, 3, 2,
     0, 0,-1,-1,-1,-1,-1,-1,-1, 3,-1, 1,-1,-1,-1,-1, 2,-1,-1,-1, 1,-1,-1,-1, 2, 3,-1, 1, 1,-1, 3, 2
};

signed char const lut5[1024] = {
     2, 2,-1, 4, 3,-1, 3, 4,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1,-1,-1,-1, 4, 3,
     1, 1,-1,-1,-1,-1,-1,-1,-1, 4,-1, 2,-1,-1,-1,-1, 3,-1,-1,-1, 2,-1,-1,-1, 3, 4,-1, 2, 2,-1, 4, 3,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1, 3,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1,-1,-1,-1, 4,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1, 2,-1, 4,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     1, 1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1, 3, 4,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4, 3,-1, 3, 4,-1, 2, 2,
     2, 2,-1, 4, 3,-1, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1, 4, 3,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1, 1, 1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 4,-1, 2,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 4,-1,-1,-1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     3,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     3, 4,-1, 2, 2,-1, 4, 3,-1,-1,-1, 2,-1,-1,-1, 3,-1,-1,-1,-1, 2,-1, 4,-1,-1,-1,-1,-1,-1,-1, 1, 1,
     3, 4,-1,-1,-1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1, 4, 3,-1, 3, 4,-1, 2, 2
};
//@}

///@internal Finds the median of 4 points along dimension @c dim.
template <int K, typename DataT>
inline int median4(Point<K, DataT> * const points, int const dim) {
    double v0 = points[0].coords.coeff(dim);
    double v1 = points[1].coords.coeff(dim);
    double v2 = points[2].coords.coeff(dim);
    double v3 = points[3].coords.coeff(dim);
    // avoids branches, but always requires 6 comparisons rather than
    // the 3-5 required by the optimal sorting network approach.
    int i = (static_cast<int>(v0 < v1) << 5) |
            (static_cast<int>(v0 < v2) << 4) |
            (static_cast<int>(v0 < v3) << 3) |
            (static_cast<int>(v1 < v2) << 2) |
            (static_cast<int>(v1 < v3) << 1) |
            static_cast<int>(v2 < v3);
    return lut4[i];
}

///@internal Finds the median of 5 points along dimension @c dim.
template <int K, typename DataT>
inline int median5(Point<K, DataT> * const points, int const dim) {
    double v0 = points[0].coords.coeff(dim);
    double v1 = points[1].coords.coeff(dim);
    double v2 = points[2].coords.coeff(dim);
    double v3 = points[3].coords.coeff(dim);
    double v4 = points[4].coords.coeff(dim);
    // avoids branches but always performs 10 comparisons, whereas
    // an optimal sorting network approach would take at most 9.
    int i = (static_cast<int>(v0 < v1) << 9) |
            (static_cast<int>(v0 < v2) << 8) |
            (static_cast<int>(v0 < v3) << 7) |
            (static_cast<int>(v0 < v4) << 6) |
            (static_cast<int>(v1 < v2) << 5) |
            (static_cast<int>(v1 < v3) << 4) |
            (static_cast<int>(v1 < v4) << 3) |
            (static_cast<int>(v2 < v3) << 2) |
            (static_cast<int>(v2 < v4) << 1) |
             static_cast<int>(v3 < v4);
    return lut5[i];
}

/** @internal
  * Returns the index of the "median of medians" for an array. This primitive
  * is used for pivot selection in the median finding algorithm.
  *
  * @pre    points != 0
  * @pre    numPoints > 0
  * @pre    dim >= 0 && dim < K
  */
template <int K, typename DataT>
int medianOfMedians(Point<K, DataT> * const points,
                    int const numPoints,
                    int const dim)
{
    int n = numPoints;
    int m = 0;
    while (true) {
        if (n <= 5) {
            switch (n) {
                case 1: break;
                case 2: m = median2(points, dim); break;
                case 3: m = median3(points, dim); break;
                case 4: m = median4(points, dim); break;
                case 5: m = median5(points, dim); break;
            }
            break;
        }
        int j = 0;
        for (int i = 0; i < n - 4; i += 5, ++j) {
            std::swap(points[j], points[i + median5(points + i, dim)]);
        }
        n = j;
    }
    return m;
}

/** @internal
  * Finds the median element of the input point array along dimension @c dim
  * and partitions the array around this element. The implementation uses the
  * linear-time "median-of-medians" algorithm.
  *
  * @pre    points != 0
  * @pre    numPoints > 0
  * @pre    dim >= 0 && dim < K
  *
  * @return The index of the median element - always equal to @c numPoints/2
  *         (rounded down).
  */
template <int K, typename DataT>
int median(Point<K, DataT> * const points,
           int const numPoints,
           int const dim)
{
    int k = numPoints >> 1;
    int left = 0;
    int right = numPoints;
    while (true) {
        int i = medianOfMedians(points + left, right - left, dim);
        i = left + partition(points + left, right - left, dim, i);
        if (k == i) {
            return k;
        } else if (k < i) {
            right = i;
        } else {
            left = i + 1;
        }
    }
}

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
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "pointer to Point array is null");
    }
    if (numPoints <= 0) {
        throw LSST_EXCEPT(except::InvalidParameterException,
                          "number of input points must be > 0");
    }
    if (pointsPerLeaf < 1) {
        throw LSST_EXCEPT(except::InvalidParameterException,
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
                right = left + median(_points + left, right - left, extDim.second);
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

}}}} // namespace lsst::ap::cluster::optics

#endif // LSST_AP_CLUSTER_OPTICS_KDTREE_CC
