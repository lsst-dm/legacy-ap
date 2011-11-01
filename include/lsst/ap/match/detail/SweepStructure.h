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
  * @brief Sweep line data structures useful for region-region matching.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#ifndef LSST_AP_MATCH_DETAIL_SWEEPSTRUCTURE_H
#define LSST_AP_MATCH_DETAIL_SWEEPSTRUCTURE_H

#include <vector>

#include "../BBox.h"
#include "../../utils/Arena.h"


namespace lsst { namespace ap { namespace match { namespace detail {

// -- Interval tree nodes ----

/** Nodes in an interval tree are either red or black.
  */
enum Color {
    BLACK = 0, RED
};

/** A node in an x-interval tree for a 2 dimensional cartesian space.
  */
struct CartesianNode {
    CartesianNode *link[2];
    Color color;
    double reach;
    double minCoord0;
    double maxCoord0;
    BBox *bbox;

    inline CartesianNode();
    inline CartesianNode(BBox *b);
    inline ~CartesianNode();
};

inline bool lessThan(CartesianNode const *left, CartesianNode const *right);

inline bool operator<(std::pair<double, CartesianNode *> const &left,
                      std::pair<double, CartesianNode *> const &right);


/** A node in a theta (longitude/right ascension) interval tree for a
  * 2 dimensional spherical space.
  */
struct SphericalNode {
    SphericalNode *link[2];
    Color color;
    unsigned int foundBy;
    double reach;
    double minCoord0;
    double maxCoord0;
    BBox *bbox;
    SphericalNode *twin;

    inline SphericalNode();
    inline SphericalNode(BBox *b, double minC0, double maxC0);
    inline ~SphericalNode();

    inline bool markFound(unsigned int searchId);
};

inline bool lessThan(SphericalNode const *left, SphericalNode const *right);

inline bool operator<(std::pair<double, SphericalNode *> const &left,
                      std::pair<double, SphericalNode *> const &right);


// -- Sweep structures ----

/** A base class for sweep-structures, templated on tree node type.
  */
template <typename Node>
class SweepStructure {
public:
    SweepStructure() : _arena(), _root(0), _heap() { }
    virtual ~SweepStructure();

    inline size_t size() const;
    inline bool empty() const;
    inline size_t getNumBytes() const;

    virtual bool isValid() const;

protected:
    typedef std::pair<double, Node *> HeapEntry;

    void _insert(Node *i);
    void _insert(BBox *b);
    void _remove(Node *n);
    void _grow(typename std::vector<HeapEntry>::size_type n);

    lsst::ap::utils::Arena<Node> _arena;
    Node *_root; ///< Root of the interval tree.
    std::vector<HeapEntry> _heap; ///< Min-heap on maximum coordinate-1 value.

private:
    // Disable copy construction and assignment
    SweepStructure(SweepStructure const &);
    SweepStructure & operator = (SweepStructure const &);

    // Low-level implementation
    static inline double _reach(Node const *n);
    inline Node * _rotate(Node *n, int dir);
    inline Node * _doubleRotate(Node *n, int dir);

    // Invariant checks
    static bool _isBinarySearchTree(Node const *n);
    static bool _isRedChildBlack(Node const *n);
    static bool _isBalanced(Node const *n, int numBlack);
    static double _checkReach(Node const *n);
    bool _isBalanced() const;
    bool _isIntervalTree() const;
};


/** Tracks a set of 2D boxes in a cartesian coordinate system
  * overlapping a sweep-line.
  *
  * @par
  * The Region type must implement the lsst::ap::match::BBox interface,
  * i.e. provide:
  *
  *   @li <tt> double getMinCoord0() const; </tt>
  *   @li <tt> double getMaxCoord0() const; </tt>
  *   @li <tt> double getMinCoord1() const; </tt>
  *   @li <tt> double getMaxCoord1() const; </tt>
  *
  * In what follows, coordinate 0 is referred to as x and coordinate 1 as y;
  * the sweep-line is a line of constant y.
  *
  * @par
  * The sweep structure supports O(log N + K) time retrieval of K boxes
  * intersecting a given x range and O(K log N) time removal of the K lowest
  * (in maximum y) boxes, where N is the number of boxes in the structure.
  * Both insertion and deletion of a single box is supported in O(log N) time.
  *
  * @par
  * The data structure used to accomplish this is an interval tree,
  * implemented by taking a red-black tree on the minimum x value
  * of member boxes and augmenting each node S with the maximum x value
  * of any box in the tree rooted at S. A separate min-heap on maximum y
  * is used to support efficient location of the node containing the box
  * with lowest maximum y.
  *
  * @par
  * This data structure is @b not thread-safe. The regions inserted
  * into the tree are @b not owned by this class and must not be
  * destroyed while they are in the sweep structure.
  */
template <typename Region>
class CartesianSweep : public SweepStructure<CartesianNode> {
public:
    CartesianSweep() : SweepStructure<CartesianNode>() { }
    virtual ~CartesianSweep() { }

    inline void insert(Region *region);

    template <typename RegionProcessor>
    void advance(double to, RegionProcessor &f);

    template <typename RegionProcessor>
    void clear(RegionProcessor &f);

    template <typename OtherRegion, typename MatchProcessor>
    void search(OtherRegion *r, MatchProcessor &f);
};


/** Tracks a set of 2D boxes in a spherical coordinate system
  * overlapping a sweep-line.
  *
  * @par
  * The Region type must implement the BBox interface, i.e. provide:
  *
  *   @li <tt> double getMinCoord0() const; </tt>
  *   @li <tt> double getMaxCoord0() const; </tt>
  *   @li <tt> double getMinCoord1() const; </tt>
  *   @li <tt> double getMaxCoord1() const; </tt>
  *
  * Coordinate 0 must correspond to a longitude/right ascension (theta) and
  * coordinate 1 to a latitude/declination (phi); the sweep-line is a line
  * of constant phi. The theta interval need not be range reduced; e.g. valid
  * theta intervals for a box 2 radians wide centered on
  * (theta, phi) = (0, 0) are:
  *
  *   @li <tt> [-1.0, 1.0] </tt>
  *   @li <tt> [2*M_PI - 1, 2*M_PI + 1] </tt>
  *   @li <tt> [2*M_PI - 1, 1] </tt>
  *
  * Coordinate values are assumed to be in units of radians.
  *
  * @par
  * The sweep structure supports O(log N + K) time retrieval of K boxes
  * intersecting a given theta range and O(K log N) time removal of the
  * K lowest (in maximum phi) boxes, where N is the number of boxes in
  * the structure. Both insertion and deletion of a single box is supported
  * in O(log N) time.
  *
  * @par
  * The data structure used to accomplish this is an interval tree,
  * implemented by taking a red-black tree on the minimum theta value
  * of member boxes and augmenting each node S with the maximum theta
  * of any box in the tree rooted at S. A separate min-heap on maximum
  * phi is used to support efficient location of the node containing the
  * box with lowest maximum phi.
  *
  * @par
  * As previously pointed out, a box B can wrap across the 0/2*M_PI radian
  * theta discontinuity. This situation is handled by inserting twin nodes
  * for such boxes: one for the theta interval <tt> [0.0, max] </tt> and
  * another for <tt> [min, 2*M_PI] </tt>.
  *
  * @par
  * Since a box may be represented by more than one node in the underlying
  * tree, and because a theta search interval can itself wrap across
  * the 0/2*M_PI radian discontinuity, care must be taken to avoid reporting
  * a box more than once during a search. This is accomplished by storing a
  * pointer to a (possibly null) "twin" node in each node, along with an
  * integer identifier for the search call that last reported the associated
  * box. Before reporting the box B associated with a node S, the search
  * identifiers for S and its twin are checked to see whether or not B
  * has already been reported by the search in progress.
  *
  * @par
  * Strictly speaking, this makes searches O(N), since the search id for
  * every node in the tree must be reset once the search id for the tree
  * (an integer) wraps. However, this is only necessary once every
  * 2^32 - 1 search calls, so in practice the logarithmic bounds hold.
  *
  * @par
  * This data structure is @b not thread-safe. The regions inserted
  * into the tree are @b not owned by this class and must not be
  * destroyed while they are in the sweep structure.
  */
template <typename Region>
class SphericalSweep : public SweepStructure<SphericalNode> {
public:
    SphericalSweep() : SweepStructure<SphericalNode>(), _searchId(1u) { }
    virtual ~SphericalSweep() { }

    inline void insert(Region *region);

    template <typename RegionProcessor>
    void advance(double to, RegionProcessor &f);

    template <typename RegionProcessor>
    void clear(RegionProcessor &f);

    template <typename OtherRegion, typename MatchProcessor>
    void search(OtherRegion *r, MatchProcessor &f);

    virtual bool isValid() const;

    void setSearchId(unsigned int searchId);

private:
    template <typename OtherRegion, typename MatchProcessor>
    void _search(OtherRegion *r, MatchProcessor &f, double min, double max);

    unsigned int _searchId; ///< Id for next search() call - must be greater
                            ///  than the search id of any tree node.
};

}}}} // namespace lsst::ap::match::detail


#include "SweepStructure.cc"

#endif // LSST_AP_MATCH_DETAIL_SWEEPSTRUCTURE_H

