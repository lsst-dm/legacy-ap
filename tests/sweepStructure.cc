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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SweepStructure

#include "boost/test/unit_test.hpp"

#include <cassert>
#include <algorithm>
#include <limits>
#include <set>
#include <vector>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/math/Random.h"

#include "lsst/ap/match/detail/SweepStructure.h"


using std::min;
using std::max;
using std::numeric_limits;
using std::set;
using std::swap;
using std::vector;

using lsst::pex::exceptions::Exception;
using lsst::afw::math::Random;
using lsst::ap::match::BBox;
using lsst::ap::match::detail::SphericalSweep;
using lsst::ap::match::detail::CartesianSweep;


/** Test region - a simple rectangle.
  */
class Box : public BBox {
public:
    Box(double min0, double max0, double min1, double max1, unsigned int id) :
        _min0(min0), _max0(max0), _min1(min1), _max1(max1), _id(id)
    { }

    virtual double getMinCoord0() const { return _min0; }
    virtual double getMaxCoord0() const { return _max0; }
    virtual double getMinCoord1() const { return _min1; }
    virtual double getMaxCoord1() const { return _max1; }
    unsigned int getId() const { return _id; }

private:
    double _min0;
    double _max0;
    double _min1;
    double _max1;
    unsigned int _id;
};

/** Do-nothing callbacks on sweep structure events.
  */
struct NoopCallbacks {
    void operator()(Box const *) const { }
    void operator()(Box const *, Box const *) { }
};

/** Callbacks for sweep structure events.
  */
class Callbacks {
public:
    bool empty() const { return _expected.empty(); }

    void expect(unsigned int i) { _expected.insert(i); }
    void clear() { _expected.clear(); }

    // called on box removal (sweep advance/clear)
    void operator()(Box const *box) {
        BOOST_CHECK_EQUAL(_expected.erase(box->getId()), 1u);
    }

    // called on match pair (sweep search)
    void operator()(Box const *box, Box const *match) {
        BOOST_CHECK_EQUAL(_expected.erase(match->getId()), 1u);
    }

private:
    set<unsigned int> _expected;
};


/** Stores a random permutation of the integers 0, 1, ..., n - 1 in v.
  */
void randomPermutation(unsigned int n, vector<unsigned int> &v, Random &rng) {
    v.clear();
    for (unsigned int i = 0; i < n; ++i) {
        v.push_back(i);
    }
    for (unsigned int i = 0; i < 2*n; ++i) {
        swap(v[rng.uniformInt(n)], v[rng.uniformInt(n)]);
    }
}

/** Determines expected matches by brute force.
  */
void findMatches(Callbacks &c, Box &b, vector<Box> &boxes) {
    for (vector<Box>::iterator i = boxes.begin(), e = boxes.end(); i != e; ++i) {
        if (i->getMinCoord0() > b.getMaxCoord0() ||
            i->getMaxCoord0() < b.getMinCoord0()) {
            continue;
        }
        c.expect(i->getId());
    }
}

/** Determines expected sweep structure removals by brute force.
  */
void findRemovals(Callbacks &c,
                  vector<Box> const &boxes,
                  double const start, ///< previous watermark
                  double const to     ///< watermark to advance to
                 )
{
    for (vector<Box>::const_iterator i = boxes.begin(), e = boxes.end();
         i != e; ++i) {
        if (i->getMaxCoord1() < start) {
            continue;
        } else if (i->getMaxCoord1() < to) {
            c.expect(i->getId());
        }
    }
}


// Explicit sweep structure instantiations
/// @cond
namespace lsst { namespace ap { namespace match { namespace detail {
    template class CartesianSweep<Box>;
    template class SphericalSweep<Box>;
}}}}
/// @endcond


template <typename Sweep>
void sweepTest() {
    static unsigned int const NBOXES = 100;

    Random rng;
    vector<unsigned int> indexes;
    vector<Box> boxes;

    // create test boxes
    for (unsigned int i = 0; i < NBOXES; ++i) {
        boxes.push_back(Box(i, i + 1.1, i, i + 1.0, i));
    }

    for (unsigned int test = 0; test < 20; ++test) {
        Sweep sweep;
        if (test == 0) {
            // insert in order of increasing coordinate-0
            for (unsigned int i = 0; i < NBOXES; ++i) {
                sweep.insert(&boxes[i]);
                BOOST_CHECK(sweep.isValid());
            }
        } else if (test == 1u) {
            // insert in order of decreasing coordinate-0
            for (unsigned int i = NBOXES - 1; i > 0; --i) {
                sweep.insert(&boxes[i]);
                BOOST_CHECK(sweep.isValid());
            }
            sweep.insert(&boxes[0]);
            BOOST_CHECK(sweep.isValid());
        } else {
            // insert in random order
            randomPermutation(NBOXES, indexes, rng);
            for (unsigned int i = 0; i < NBOXES; ++i) {
                sweep.insert(&boxes[indexes[i]]);
                BOOST_CHECK(sweep.isValid());
            }
        }
        Callbacks c;
        for (unsigned int i = 0; i < NBOXES; ++i) {
            c.expect(i);
            if (i > 0) {
                c.expect(i - 1);
            }
            if (i < NBOXES - 1) {
                c.expect(i + 1);
            }
            sweep.search(&boxes[i], c);
            BOOST_CHECK(c.empty());
            c.clear();
        }
        for (unsigned int i = 0; i < NBOXES; ++i) {
            c.expect(i);
            sweep.advance(boxes[i].getMaxCoord1() + 0.1, c);
            BOOST_CHECK(sweep.isValid());
            BOOST_CHECK(c.empty());
            c.clear();
        }
        BOOST_CHECK_EQUAL(sweep.size(), 0u);
    }
}

BOOST_AUTO_TEST_CASE(cartesianSweep) {
    sweepTest<CartesianSweep<Box> >();
}
BOOST_AUTO_TEST_CASE(sphericalSweep) {
    sweepTest<SphericalSweep<Box> >();
}

template <typename Sweep>
void sweepRandomTest() {
    static unsigned int const NBOXES = 20;

    Random rng;
    vector<Box> boxes;
    Sweep sweep;
    Callbacks c;
    NoopCallbacks nc;

    // create test boxes
    for (unsigned int test = 0; test < 100; ++test) {
        for (unsigned int i = 0; i < NBOXES; ++i) {
            double c00 = rng.uniform();
            double c01 = rng.uniform();
            double c10 = rng.uniform();
            double c11 = rng.uniform();
            boxes.push_back(Box(min(c00, c01), max(c00, c01),
                                min(c10, c11), max(c10, c11), i));
        }
        for (unsigned int i = 0; i < NBOXES; ++i) {
            sweep.insert(&boxes[i]);
            BOOST_CHECK(sweep.isValid());
        }
        for (unsigned int i = 0; i < NBOXES; ++i) {
            findMatches(c, boxes[i], boxes);
            sweep.search(&boxes[i], c);
            BOOST_CHECK(c.empty());
            c.clear();
        }
        double start = -numeric_limits<double>::infinity();
        for (unsigned int i = 0; i < NBOXES; ++i) {
            double to = max(start, boxes[i].getMaxCoord1() + 0.1);
            findRemovals(c, boxes, start, to);
            sweep.advance(to, c);
            BOOST_CHECK(sweep.isValid());
            BOOST_CHECK(c.empty());
            c.clear();
            start = to;
        }
        BOOST_CHECK_EQUAL(sweep.size(), 0u);
        boxes.clear();
        sweep.clear(nc);
    }
}

BOOST_AUTO_TEST_CASE(cartesianSweepRandom) {
    sweepRandomTest<CartesianSweep<Box> >();
}
BOOST_AUTO_TEST_CASE(sphericalSweepRandom) {
    sweepRandomTest<SphericalSweep<Box> >();
}

unsigned int const FACTORIAL[] = {
    0,
    1,
    2,
    2*3,
    2*3*4,
    2*3*4*5,
    2*3*4*5*6,
    2*3*4*5*6*7
};

void permutation(unsigned int n, unsigned int pi, unsigned int pr, vector<Box> &v) {
    assert(n > 0);
    v.clear();
    vector<unsigned int> inserts;
    vector<unsigned int> removes;
    for (unsigned int i = 0; i < n; ++i) {
        inserts.push_back(i);
        removes.push_back(i);
    }
    for (unsigned int i = n; i > 0; --i) {
        unsigned int which = pi % i;
        pi /= i;
        unsigned int c0 = inserts[which];
        inserts.erase(inserts.begin() + which);
        which = pr % i;
        pr /= i;
        unsigned int c1 = removes[which];
        removes.erase(removes.begin() + which);
        v.push_back(Box(c0, c0 + 1, c1, c1 + 1, i));
    }
}

template <typename Sweep>
void sweepCornerCasesTest(unsigned int maxNodes) {
    Callbacks c;
    NoopCallbacks nc;
    Sweep sweep;
    Box b1(0.0, 0.0, 0.0, 0.0, 0u);
    Box b2(1.0, 1.0, 0.0, 0.0, 1u);
    vector<Box> boxes;
    maxNodes = min(maxNodes, static_cast<unsigned int>(
                   sizeof(FACTORIAL)/sizeof(unsigned int) - 1));

    // 1. test search on empty tree
    sweep.search(&b1, c);
    // 2. test advance on empty tree
    sweep.advance(1.0, c); 
    // 3. test clear on empty tree
    sweep.clear(c);
    // 4. test search that should not return anything
    sweep.insert(&b1);
    sweep.search(&b2, c);
    // 5. test all possible tree arrangements with 1-6 nodes
    for (unsigned int n = 1; n <= maxNodes; ++n) {
        for (unsigned int pi = 0; pi < FACTORIAL[n]; ++pi) {
            for (unsigned int pr = 0; pr < FACTORIAL[n]; ++pr) {
                sweep.clear(nc);
                permutation(n, pi, pr, boxes);
                for (unsigned int i = 0; i < n; ++i) {
                    sweep.insert(&boxes[i]);
                    BOOST_CHECK(sweep.isValid());
                }
                for (unsigned int i = 0; i < n; ++i) {
                    sweep.advance(i + 1.1, nc);
                    BOOST_CHECK(sweep.isValid());
                }
                BOOST_CHECK(sweep.empty());
            }
        }
    }
    // 6. test trees containing nodes with identical spatial extents 
    for (unsigned int n = 2; n <= 10; ++n) {
        sweep.clear(nc);
        boxes.clear();
        for (unsigned int i = 0; i < n; ++i) {
            boxes.push_back(Box(0.0, 1.0, i, i + 1, i));
        }
        for (unsigned int i = 0; i < n; ++i) {
            sweep.insert(&boxes[i]);
            BOOST_CHECK(sweep.isValid());
        }
        for (unsigned int i = 0; i < n; ++i) {
            sweep.advance(i + 1.1, nc);
            BOOST_CHECK(sweep.isValid());
        }
        BOOST_CHECK(sweep.empty());
    }
}

BOOST_AUTO_TEST_CASE(cartesianSweepCornerCases) {
    sweepCornerCasesTest<CartesianSweep<Box> >(5u);
}
BOOST_AUTO_TEST_CASE(sphericalSweepCornerCases) {
    sweepCornerCasesTest<SphericalSweep<Box> >(5u);
}


// Sweep structures over the sphere have additional subtleties

BOOST_AUTO_TEST_CASE(sphericalWrap) {
    Box search1(0.0, 20.0, 0.0, 0.0, 0);
    Box search2(340.0, 359.9999, 0.0, 0.0, 0);
    Box searchWrap(-20.0, 20.0, 0.0, 0.0, 0);
    Box b1(350.0, 10.0, 0.0, 1.0, 0);
    Box b2(350.0, 370.0, 0.0, 1.0, 1);
    Box b3(-10.0, 10.0, 0.0, 1.0, 2);
    Box b4(0.0, 5.0, 1.0, 2.0, 3);
    Box b5(355.0, 359.9999, 2.0, 3.0, 4);
    Box b6(180.0, 181.0, 2.0, 3.0, 5);

    // Test wrapping regions
    Callbacks c;
    NoopCallbacks nc;
    SphericalSweep<Box> sweep;
    sweep.insert(&b1); BOOST_CHECK(sweep.isValid());
    sweep.insert(&b2); BOOST_CHECK(sweep.isValid());
    sweep.insert(&b3); BOOST_CHECK(sweep.isValid());
    sweep.insert(&b4); BOOST_CHECK(sweep.isValid());
    sweep.insert(&b5); BOOST_CHECK(sweep.isValid());
    sweep.insert(&b6); BOOST_CHECK(sweep.isValid());
    c.expect(0); c.expect(1); c.expect(2); c.expect(3);
    // ... with search intervals that do not wrap
    sweep.search(&search1, c);
    BOOST_CHECK(c.empty());
    c.clear();
    c.expect(0); c.expect(1); c.expect(2); c.expect(4);
    sweep.search(&search2, c);
    BOOST_CHECK(c.empty());
    // ... with a search interval that also wraps
    c.clear();
    c.expect(0); c.expect(1); c.expect(2); c.expect(3); c.expect(4);
    sweep.search(&searchWrap, c);
    BOOST_CHECK(c.empty());
    c.clear();
    c.expect(0); c.expect(1); c.expect(2);
    sweep.advance(1.1, c);
    BOOST_CHECK(sweep.isValid());
    BOOST_CHECK(c.empty());
    c.clear();
    c.expect(3); c.expect(4); c.expect(5);
    sweep.clear(c);
    BOOST_CHECK(sweep.isValid());
    BOOST_CHECK(sweep.empty());
    BOOST_CHECK(c.empty());

    // Test wrapping search intervals on trees without any wrapping regions
    sweep.clear(nc);
    c.clear();
    sweep.insert(&b4); BOOST_CHECK(sweep.isValid());
    sweep.insert(&b5); BOOST_CHECK(sweep.isValid());
    c.expect(3); c.expect(4);
    sweep.search(&searchWrap, c);
    BOOST_CHECK(c.empty());
    c.expect(3); c.expect(4);
    sweep.advance(3.1, c);
    BOOST_CHECK(sweep.isValid());
    BOOST_CHECK(sweep.empty());
    BOOST_CHECK(c.empty());
}

BOOST_AUTO_TEST_CASE(sphericalSearchIdWrap) {
    // Test that search id wrapping is handled correctly
    Box b(10.0, 20.0, 0.0, 0.0, 0);
    Callbacks c;
    SphericalSweep<Box> sweep;
    sweep.insert(&b);
    // make search id wrap to 0 on next call to search
    sweep.setSearchId(numeric_limits<unsigned int>::max());
    BOOST_CHECK(sweep.isValid());
    c.expect(0);
    sweep.search(&b, c);
    BOOST_CHECK(sweep.isValid());
    BOOST_CHECK(c.empty());
    c.expect(0);
    sweep.search(&b, c);
    BOOST_CHECK(c.empty());
}

BOOST_AUTO_TEST_CASE(nodeSize) {
    // not a test-case - this just prints out the size of a
    // cartesian/spherical tree node on the target platform.
    BOOST_MESSAGE("sizeof(CartesianNode) = " <<
                  sizeof(lsst::ap::match::detail::CartesianNode));
    BOOST_MESSAGE("sizeof(SphericalNode) = " <<
                  sizeof(lsst::ap::match::detail::SphericalNode));
}

