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
  * @brief Inline/templated sweep-line data structure method implementations.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_MATCH_DETAIL_SWEEPSTRUCTURE_CC
#define LSST_AP_MATCH_DETAIL_SWEEPSTRUCTURE_CC

#include "SweepStructure.h"

#include <math.h>
#include <climits>
#include <algorithm>

#include "lsst/afw/geom/Angle.h"
#include "lsst/pex/exceptions.h"
#include "../../utils/SpatialUtils.h"

namespace afwGeom = lsst::afw::geom;

namespace lsst { namespace ap { namespace match { namespace detail {

// -- CartesianNode implementation ----

/** Creates a dummy tree root.
  */
inline CartesianNode::CartesianNode() :
    color(BLACK),
    reach(0.0),
    minCoord0(0.0),
    maxCoord0(0.0),
    bbox(0)
{
    link[0] = 0;
    link[1] = 0;
}

/** Creates a new CartesianNode from the given BBox.
  * Ownership of the given BBox is @b not transferred.
  */
inline CartesianNode::CartesianNode(BBox *b) :
    color(RED),
    reach(b->getMaxCoord0()),
    minCoord0(b->getMinCoord0()),
    maxCoord0(reach),
    bbox(b)
{
    link[0] = 0;
    link[1] = 0;
}

inline CartesianNode::~CartesianNode() {
#ifndef NDEBUG
    link[0] = 0;
    link[1] = 0;
    bbox = 0;
#endif
}

/** Orders tree nodes by minimum coordinate-0 (x) value, using
  * the addresses of embedded BBox instances to break ties.
  */
inline bool lessThan(CartesianNode const *a, CartesianNode const *b) {
    return (a->minCoord0 == b->minCoord0) ?
           a->bbox < b->bbox : a->minCoord0 < b->minCoord0;
}

inline bool operator<(std::pair<double, CartesianNode *> const &n1,
                      std::pair<double, CartesianNode *> const &n2)
{
    return n1.first > n2.first;
}


// -- SphericalNode implementation ----

/** Creates a dummy tree root.
  */
inline SphericalNode::SphericalNode() :
    color(BLACK),
    foundBy(0),
    reach(0.0),
    minCoord0(0.0),
    maxCoord0(0.0),
    bbox(0),
    twin(0)
{
    link[0] = 0;
    link[1] = 0;
}

/** Creates a new SphericalNode from the given BBox and theta interval.
  * Ownership of the given BBox is @b not transferred.
  */
inline SphericalNode::SphericalNode(BBox *b,
                                    double minC0,
                                    double maxC0
                                   ) :
    color(RED),
    foundBy(0),
    reach(maxC0),
    minCoord0(minC0),
    maxCoord0(maxC0),
    bbox(b),
    twin(0)
{
    link[0] = 0;
    link[1] = 0;
}

inline SphericalNode::~SphericalNode() {
    if (twin != 0) {
        twin->twin = 0;
    }
#ifndef NDEBUG
    link[0] = 0;
    link[1] = 0;
    bbox = 0;
    twin = 0;
#endif
}

inline bool SphericalNode::markFound(unsigned int searchId) {
    unsigned int id = foundBy;
    foundBy = searchId;
    return (id == searchId || (twin != 0 && twin->foundBy == searchId));
}

/** Orders tree nodes by minimum coordinate-0 (longitude/right-ascension)
  * value, using the addresses of embedded BBox instances to break ties.
  */
inline bool lessThan(SphericalNode const *a, SphericalNode const *b) {
    return (a->minCoord0 == b->minCoord0) ?
           a->bbox < b->bbox : a->minCoord0 < b->minCoord0;
}

inline bool operator<(std::pair<double, SphericalNode *> const &left,
                      std::pair<double, SphericalNode *> const &right)
{
    return left.first > right.first;
}


// -- SweepStructure implementation ---

/** Returns the number of nodes in the sweep structure.
  */
template <typename Node>
inline size_t SweepStructure<Node>::size() const {
    return _heap.size();
}

/** Returns true if the sweep structure is empty.
  */
template <typename Node>
inline bool SweepStructure<Node>::empty() const {
    return size() == 0;
}

/** Returns an estimate for the amount of memory (in bytes) allocated
  * by the sweep-structure.
  */
template <typename Node>
inline size_t SweepStructure<Node>::getNumBytes() const {
    return _arena.getNumBytes() + _heap.capacity()*sizeof(HeapEntry);
}

/** Finds the maximum coordinate-0 value of any node in the tree rooted at n
  * (its reach) by comparing the maximum coordinate-0 value of the
  * coordinate-0 interval associated with n to the reach of both children.
  */
template <typename Node>
inline double SweepStructure<Node>::_reach(Node const *n) {
    double reach = n->maxCoord0;
    Node const *c0 = n->link[0];
    if (c0 != 0 && c0->reach > reach) {
        reach = c0->reach;
    }
    Node const *c1 = n->link[1];
    if (c1 != 0 && c1->reach > reach) {
        reach = c1->reach;
    }
    return reach;
}

/** Performs a tree rotation.
  */
template <typename Node>
inline Node * SweepStructure<Node>::_rotate(Node *n, int dir) {
    Node *r = n->link[!dir];
    n->link[!dir] = r->link[dir];
    r->link[dir] = n;
    r->color = BLACK;
    r->reach = n->reach;
    n->color = RED;
    n->reach = _reach(n);
    return r;
}

/** Performs a double tree rotation.
  */
template <typename Node>
inline Node * SweepStructure<Node>::_doubleRotate(Node *n, int dir) {
    n->link[!dir] = _rotate(n->link[!dir], !dir);
    return _rotate(n, dir); 
}

// -- CartesianSweep implementation ---

/** Inserts the given Region into the sweep structure - the sweep
  * structure does @b not assume ownership of the region.
  */
template <typename Region>
inline void CartesianSweep<Region>::insert(Region *region) {
    _insert(region);
}

/** Removes all regions R with maximum y value less than @a to from the sweep
  * structure. Invokes f(R) exactly once for each region R that is removed.
  *
  * @tparam RegionProcessor
  *         A functor providing the following operator:
  *         @code void RegionProcessor::operator()(Region &) @endcode
  *         A RegionProcessor may have state.
  */
template <typename Region>
    template <typename RegionProcessor>
void CartesianSweep<Region>::advance(double to,
                                     RegionProcessor &f
                                    )
{
    while (!_heap.empty() && _heap.front().first < to) {
        std::pop_heap(_heap.begin(), _heap.end());
        CartesianNode *n = _heap.back().second;
        _heap.pop_back();
        _remove(n);
        f(dynamic_cast<Region *>(n->bbox));
        _arena.destroy(n);
    }
}

/** Empties the sweep structure, invoking f(R) exactly once for each
  * region R that is removed.
  *
  * @tparam RegionProcessor
  *         A functor providing the following operator:
  *         @code void RegionProcessor::operator()(Region &) @endcode
  *         A RegionProcessor may have state.
  */
template <typename Region>
    template <typename RegionProcessor>
void CartesianSweep<Region>::clear(RegionProcessor &f) {
    typedef std::vector<HeapEntry>::iterator Iter;
    for (Iter i = _heap.begin(), e = _heap.end(); i != e; ++i) {
        CartesianNode *n = i->second;
        f(dynamic_cast<Region *>(n->bbox));
        _arena.destroy(n);
    }
    _heap.clear();
    _root = 0;
}

/** Invokes f(r, R) exactly once for each region R
  * in the sweep structure that overlaps the x-interval of r.
  *
  * @tparam OtherRegion
  *         A class that provides the lsst::ap::match::BBox interface,
  *         though inheritance from the interface is not required.
  * @tparam MatchProcessor
  *         A functor providing the following operator:
  *         @code void MatchProcessor::operator()(OtherRegion &, Region &) @endcode
  *         A MatchProcessor may have state.
  */
template <typename Region>
    template <typename OtherRegion, typename MatchProcessor>
void CartesianSweep<Region>::search(OtherRegion *r,
                                    MatchProcessor &f
                                   )
{
    CartesianNode *node[2*sizeof(size_t)*CHAR_BIT];
    int dir[2*sizeof(size_t)*CHAR_BIT];
    CartesianNode *n = _root;
    if (n == 0) {
        return;
    }
    double const min = r->getMinCoord0();
    double const max = r->getMaxCoord0();
    int s = 0;
    
    while (true) {
        if (n->reach >= min) {
            if (n->minCoord0 <= max && n->maxCoord0 >= min) {
                f(r, dynamic_cast<Region *>(n->bbox));
            }
            if (n->link[0] != 0) {
                node[s] = n;
                dir[s] = 0;
                ++s;
                n = n->link[0];
                continue;
            } else if (n->link[1] != 0) {
                node[s] = n;
                dir[s] = 1;
                ++s;
                n = n->link[1];
                continue;
            }
        }
        for (--s; s >= 0 && (dir[s] == 1 || node[s]->link[1] == 0); --s) { }
        if (s < 0 || node[s]->minCoord0 > max) {
            break;
        }
        dir[s] = 1;
        n = node[s]->link[1];
        ++s;
    }
}


// -- SphericalSweep implementation ----

/** Inserts the given Region into the sweep structure - the sweep
  * structure does @b not assume ownership of the region.
  */
template <typename Region>
inline void SphericalSweep<Region>::insert(Region *region) {
    _insert(region);
}

/** Removes all regions R with maximum phi (latitude/declination) less than
  * to from the sweep structure. Invokes f(R) exactly once for each region R
  * that is removed.
  *
  * @tparam RegionProcessor
  *         A functor providing the following operator:
  *         @code void RegionProcessor::operator()(Region &) @endcode
  *         A RegionProcessor may have state.
  */
template <typename Region>
    template <typename RegionProcessor>
void SphericalSweep<Region>::advance(double to,
                                     RegionProcessor &f
                                    )
{
    while (!_heap.empty() && _heap.front().first < to) {
        std::pop_heap(_heap.begin(), _heap.end());
        SphericalNode *n = _heap.back().second;
        _heap.pop_back();
        _remove(n);
        if (n->twin == 0) {
            // Note that if a node n is deleted by this loop,
            // so is its twin. Only invoke the region processor
            // on nodes that never had or no longer have a twin -
            // this avoids double invoking on the same region.
            f(dynamic_cast<Region *>(n->bbox));
        }
        _arena.destroy(n);
    }
}

/** Empties the sweep structure, invoking f(R) exactly once for each
  * region R that is removed.
  *
  * @tparam RegionProcessor
  *         A functor providing the following operator:
  *         @code void RegionProcessor::operator()(Region &) @endcode
  *         A RegionProcessor may have state.
  */
template <typename Region>
    template <typename RegionProcessor>
void SphericalSweep<Region>::clear(RegionProcessor &f) {
    typedef std::vector<HeapEntry>::iterator Iter;
    for (Iter i = _heap.begin(), e = _heap.end(); i != e; ++i) {
        SphericalNode *n = i->second;
        if (n->twin == 0) {
            // Note that if a node n is deleted by this loop,
            // so is its twin. Only invoke the region processor
            // on nodes that never had or no longer have a twin -
            // this avoids double invoking on the same region.
            f(dynamic_cast<Region *>(n->bbox));
        }
        _arena.destroy(n);
    }
    _heap.clear();
    _root = 0;
    _searchId = 1;
}

/** Invokes f(r, R) exactly once for each region R in the sweep structure
  * that overlaps the theta (longitude/right-ascension) interval of r.
  *
  * @tparam OtherRegion
  *         A class that provides the lsst::ap::match::BBox interface,
  *         though inheritance from the interface is not required.
  * @tparam MatchProcessor
  *         A functor providing the following operator:
  *         @code void MatchProcessor::operator()(OtherRegion &, Region &) @endcode
  *         A MatchProcessor may have state.
  */
template <typename Region>
    template <typename OtherRegion, typename MatchProcessor>
void SphericalSweep<Region>::search(OtherRegion *r,
                                    MatchProcessor &f
                                   )
{
    if (_root == 0) {
        return;
    }
    double min = r->getMinCoord0();
    double max = r->getMaxCoord0();
    lsst::ap::utils::thetaRangeReduce(min, max);
    if (min > max) {
        _search(r, f, 0.0, max);
        max = afwGeom::TWOPI;
    }
    _search(r, f, min, max);
    ++_searchId;
    if (_searchId == 0) {
        // _searchId wrapped to zero. Without further action, a region
        // last reported by search N will never be reported once the search
        // id reaches N again. Therefore, whenever _searchId wraps to
        // zero, set the search id of all nodes to 0 and set _searchId to
        // 1.
        typedef std::vector<HeapEntry>::iterator Iter;
        for (Iter i = _heap.begin(), e = _heap.end(); i != e; ++i) {
            i->second->foundBy = 0;
        }
        _searchId = 1;
    }
}

template <typename Region>
bool SphericalSweep<Region>::isValid() const {
    if (!SweepStructure<SphericalNode>::isValid()) {
        return false;
    }
    // check that the search id of every node is less than the id
    // of the next search call
    typedef std::vector<HeapEntry>::const_iterator Iter;
    for (Iter i = _heap.begin(), e = _heap.end(); i != e; ++i) {
        if (i->second->foundBy >= _searchId) {
            return false;
        }
    }
    return true;
}

/** Sets the search id for the next search() call. Calling this is safe but
  * useless outside of unit-tests: search() automatically increments
  * the search id for the sweep structure.
  */
template <typename Region>
void SphericalSweep<Region>::setSearchId(unsigned int searchId) {
    typedef std::vector<HeapEntry>::iterator Iter;
    if (searchId == 0) {
        // a search id of 0 is reserved to mean "not found by any search"
        searchId = 1;
    }
    if (searchId < _searchId) {
        // the tree may contain nodes with search id equal to searchId
        for (Iter i = _heap.begin(), e = _heap.end(); i != e; ++i) {
            i->second->foundBy = 0;
        }
    }
    _searchId = searchId;
}

/** Invokes f(r, R) exactly once for each region R in the sweep structure
  * that overlaps the theta (longitude/right-ascension) interval [min, max].
  */
template <typename Region>
    template <typename OtherRegion, typename MatchProcessor>
void SphericalSweep<Region>::_search(OtherRegion *r,
                                     MatchProcessor &f,
                                     double min,
                                     double max
                                    )
{
    SphericalNode *node[2*sizeof(size_t)*CHAR_BIT];
    int dir[2*sizeof(size_t)*CHAR_BIT];
    int s = 0;
    SphericalNode *n = _root;

    while (true) {
        if (n->reach >= min) {
            if (n->minCoord0 <= max && n->maxCoord0 >= min) {
                if (!n->markFound(_searchId)) {
                    // found a match to an as yet unreported region -
                    // invoke the match processor on it
                    f(r, dynamic_cast<Region *>(n->bbox));
                }
            }
            if (n->link[0] != 0) {
                node[s] = n;
                dir[s] = 0;
                ++s;
                n = n->link[0];
                continue;
            } else if (n->link[1] != 0) {
                node[s] = n;
                dir[s] = 1;
                ++s;
                n = n->link[1];
                continue;
            }   
        }
        for (--s; s >= 0 && (dir[s] == 1 || node[s]->link[1] == 0); --s) { }
        if (s < 0 || node[s]->minCoord0 > max) {
            break;
        }
        dir[s] = 1;
        n = node[s]->link[1];
        ++s;
    }
}


}}}} // namespace lsst::ap::match::detail

#endif // LSST_AP_MATCH_DETAIL_SWEEPSTRUCTURE_CC

