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
  * @brief Sweep-line base class implementation and instantiations.
  *
  * @ingroup ap
  * @author Serge Monkewitz
  */
#include "lsst/ap/match/detail/SweepStructure.h"

#include "lsst/utils/ieee.h"


namespace lsst { namespace ap { namespace match { namespace detail {

// -- SweepStructure implementation ---

template <typename Node>
SweepStructure<Node>::~SweepStructure() {
#ifndef NDEBUG
    typedef typename std::vector<HeapEntry>::iterator Iter;
    _root = 0;
    for (Iter i = _heap.begin(), e = _heap.end(); i != e; ++i) {
        i->second = 0;
    }
#endif
}

/** Checks whether or not the data structure is valid, i.e. if
  * its theoretical invariants hold in practice. To be used for
  * unit-testing.
  */
template <typename Node>
bool SweepStructure<Node>::isValid() const {
    return _isBinarySearchTree(_root) &&
           _isRedChildBlack(_root) &&
           _isBalanced() &&
           _isIntervalTree();
}

/** Inserts the node i into the interval tree.
  */
template <typename Node>
void SweepStructure<Node>::_insert(Node *i) {
    if (_root == NULL) {
        // tree is empty
        _root = i;
        i->color = BLACK;
        return;
    }
    // top-down tree insertion algorithm
    Node head;         // dummy root
    Node *ggp = &head; // great grand parent
    Node *gp = 0;      // grand parent
    Node *p = 0;       // parent
    Node *n = _root;   // iterator
    double const maxCoord0 = i->maxCoord0;
    int dir = 0;
    int lastDir = 0;
    head.link[1] = n;

    while (true) {
        if (n == 0) {
            // insert new node
            p->link[dir] = i;
            if (p != 0 && p->color != BLACK) {
               // p and i are red: fix red violation
               int d = (ggp->link[1] == gp);
               if (i == p->link[lastDir]) {
                   ggp->link[d] = _rotate(gp, !lastDir);
               } else {
                   ggp->link[d] = _doubleRotate(gp, !lastDir);
               }
            }
            // done
            break;
        }
        // update reach
        if (maxCoord0 > n->reach) {
            n->reach = maxCoord0;
        }
        Node *c0 = n->link[0];
        Node *c1 = n->link[1];
        if (c0 != 0 && c1 != 0 && c0->color != BLACK && c1->color != BLACK) {
            // color flip
            c0->color = BLACK;
            c1->color = BLACK;
            n->color = RED;
        }
        if (n->color != BLACK && p != 0 && p->color != BLACK) {
            // fix red violation
            int d = (ggp->link[1] == gp);
            if (n == p->link[lastDir]) {
                ggp->link[d] = _rotate(gp, !lastDir);
            } else {
                ggp->link[d] = _doubleRotate(gp, !lastDir);
            }
        }
        lastDir = dir;
        dir = lessThan(n, i);
        if (gp != 0) {
            ggp = gp;
        }
        gp = p;
        p = n;
        n = n->link[dir];
    }
    _root = head.link[1];
    _root->color = BLACK;
}

/** Removes node n from the sweep structure.
  */
template <typename Node>
void SweepStructure<Node>::_remove(Node *r) {
    // top-down tree removal - push a red node down
    // the tree such that the node to be removed
    // is red. On the way down the tree, build a 
    // root to leaf path, and finally, fix up reaches
    // bottom-up.
    Node *stack[2*sizeof(size_t)*CHAR_BIT];
    Node head;       // dummy root
    Node *gp = 0;    // grand parent
    Node *p = 0;     // parent
    Node *n = &head; // iterator
    head.link[1] = _root;
    int dir = 1;
    int level = 0;

    while (n->link[dir] != 0) {
        int lastDir = dir;
        if (p != 0) {
            stack[level++] = p;
        }
        gp = p;
        p = n;
        n = n->link[dir];
        dir = lessThan(n, r);

        // push down red node
        if (n->color == BLACK &&
            (n->link[dir] == 0 || n->link[dir]->color == BLACK)) {
            if (n->link[!dir] != 0 && n->link[!dir]->color != BLACK) {
                stack[level++] = p;
                p = p->link[lastDir] = _rotate(n, dir);
            } else {
                Node *s = p->link[!lastDir]; // sibling
                if (s != 0) {
                    Node *sld = s->link[lastDir];
                    Node *snld = s->link[!lastDir];
                    if ((sld == 0 || sld->color == BLACK) &&
                        (snld == 0 || snld->color == BLACK)) {
                        // color flip
                        p->color = BLACK;
                        s->color = RED;
                        n->color = RED; 
                    } else {
                        // double or single rotation
                        int d = (gp->link[1] == p);
                        if (sld != 0 && sld->color == RED) {
                            stack[level++] = sld;
                            s = _doubleRotate(p, lastDir);
                        } else {
                            stack[level++] = s;
                            s = _rotate(p, lastDir);
                        }
                        // fixup colors
                        gp->link[d] = s;
                        n->color = RED;
                        s->color = RED;
                        s->link[0]->color = BLACK;
                        s->link[1]->color = BLACK;
                    }
                }
            }
        }
    }

    p->link[p->link[1] == n] = n->link[n->link[0] == 0];

    // walk up the node stack, recomputing reaches. Note that
    // either level == 0, or stack[0] == &head
    for (; level > 0; p = stack[--level]) {
        if (p == r) {
            // reached node to remove, exchange it with n
            n->link[0] = r->link[0];
            n->link[1] = r->link[1];
            n->color = r->color;
            gp = stack[level - 1];
            gp->link[gp->link[1] == r] = n; 
            p = n;
        }
        p->reach = _reach(p);
    }

    // update root and color it black
    _root = head.link[1];
    if (_root != 0) {
        _root->color = BLACK;
    }
}

/** Checks that a depth first traversal of the tree rooted at n yields
  * an ordered list of boxes.
  */
template <typename Node>
bool SweepStructure<Node>::_isBinarySearchTree(Node const *n) {
    if (n == 0) {
        return true;
    }
    if ((n->link[0] != 0 && lessThan(n, n->link[0])) ||
        (n->link[1] != 0 && lessThan(n->link[1], n))) {
        return false;
    }
    return _isBinarySearchTree(n->link[0]) &&
           _isBinarySearchTree(n->link[1]);
}

/** Checks that every child of a red node is black.
  */
template <typename Node>
bool SweepStructure<Node>::_isRedChildBlack(Node const *n) {
    if (n == 0) {
        return true;
    }
    if (n->color == RED) {
        if ((n->link[0] != 0 && n->link[0]->color == RED) ||
            (n->link[1] != 0 && n->link[1]->color == RED)) {
            return false;
        }
    }
    return _isRedChildBlack(n->link[0]) &&
           _isRedChildBlack(n->link[1]);
}

/** Checks that all root-to-leaf paths have the same number of black nodes.
  */
template <typename Node>
bool SweepStructure<Node>::_isBalanced() const {
    // count number of black nodes on path from root to minimum leaf
    int numBlack = 0;
    Node *n = _root;
    while (n != 0) {
        if (n->color == BLACK) {
            ++numBlack;
        }
        n = n->link[0];
    }
    // check that all other paths agree
    return _isBalanced(_root, numBlack);
}

/** Checks that every n-to-leaf path has the given number of black nodes.
  */
template <typename Node>
bool SweepStructure<Node>::_isBalanced(Node const *n, int numBlack) {
    if (n == 0) {
        return (numBlack == 0);
    }
    if (n->color == BLACK) {
        --numBlack;
    }
    return _isBalanced(n->link[0], numBlack) && 
           _isBalanced(n->link[1], numBlack);
}

/** Tests if this tree is an interval tree.
  */
template <typename Node>
bool SweepStructure<Node>::_isIntervalTree() const {
    double reach = _checkReach(_root);
    return !lsst::utils::isnan(reach);
}

/** Returns the reach of node n. If the reach of n is not equal to the
  * maximum coordinate-0 value of any node in the tree rooted at n, then
  * a quiet NaN is returned.
  */
template <typename Node>
double SweepStructure<Node>::_checkReach(Node const *n) {
    if (n == 0) {
        return -std::numeric_limits<double>::infinity();
    }
    double r0 = _checkReach(n->link[0]);
    if (lsst::utils::isnan(r0)) {
        return r0;
    }
    double r1 = _checkReach(n->link[1]);
    if (lsst::utils::isnan(r1)) {
        return r1;
    }
    double r = std::max(std::max(r0, r1), n->maxCoord0);;
    return (r == n->reach) ? r : std::numeric_limits<double>::quiet_NaN();
}


// -- SweepStructure<CartesianNode> method specializations ----

/** Inserts the given BBox into the sweep structure - note that ownership
  * of the box is @b not transferred to the sweep structure.
  */
template <>
void SweepStructure<CartesianNode>::_insert(BBox *b) {
    if (b != 0) {
        _heap.reserve(_heap.size() + 1);
        CartesianNode *n = new (_arena) CartesianNode(b);
        // after this point, no exception can be thrown
        _insert(n);
        _heap.push_back(HeapEntry(b->getMaxCoord1(), n));
        std::push_heap(_heap.begin(), _heap.end());
    }
}


// -- SweepStructure<SphericalNode> method specializations ----

/** Inserts the given BBox into the sweep structure - note that ownership
  * of the box is @b not transferred to the sweep structure.
  */
template <>
void SweepStructure<SphericalNode>::_insert(BBox *b) {
    if (b == 0) {
        return;
    }
    double min = b->getMinCoord0();
    double max = b->getMaxCoord0();
    lsst::ap::utils::thetaRangeReduce(min, max);
    if (min > max) {
        // box wraps across the 0/2*M_PI longitude angle
        // discontinuity - insert twin nodes for [0, min] and
        // [max, 2*M_PI]
        _heap.reserve(_heap.size() + 2);
        SphericalNode *n1 = new (_arena) SphericalNode(b, 0.0, max);
        SphericalNode *n2 = 0;
        try {
            n2 = new (_arena) SphericalNode(b, min, 2.0*M_PI);
        } catch (...) {
            // exception safety: if allocation for n2 fails, delete n1.
            _arena.destroy(n1);
            throw;
        }
        // after this point no exception can be thrown
        n1->twin = n2;
        n2->twin = n1;
        _insert(n1);
        _insert(n2);
        double maxCoord1 = b->getMaxCoord1();
        _heap.push_back(HeapEntry(maxCoord1, n1));
        std::push_heap(_heap.begin(), _heap.end());
        _heap.push_back(HeapEntry(maxCoord1, n2));
        std::push_heap(_heap.begin(), _heap.end());
    } else {
        // box does not wrap - insert a single node.
        _heap.reserve(_heap.size() + 1);
        SphericalNode *n = new (_arena) SphericalNode(b, min, max);
        // after this point no exception can be thrown
        _insert(n);
        _heap.push_back(HeapEntry(b->getMaxCoord1(), n));
        std::push_heap(_heap.begin(), _heap.end());
    }
}


// -- SweepStructure instantiations ---

/// @cond
template class SweepStructure<CartesianNode>;
template class SweepStructure<SphericalNode>;
/// @endcond

}}}} // namespace lsst::ap::match::detail

