/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2012 LSST Corporation.
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
 
#include <limits>

#include "boost/shared_array.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SeedList
#include "boost/test/unit_test.hpp"

#include "lsst/afw/math/Random.h"
#include "lsst/ap/cluster/optics/SeedList.cc"


namespace optics = lsst::ap::cluster::optics;

using lsst::afw::math::Random;

typedef optics::Point<1, void *> Point;
typedef optics::SeedList<1, void *> SeedList;

namespace {

boost::shared_array<Point> const makePoints(int n, bool random) {
    static Random rng(Random::MT19937);

    boost::shared_array<Point> points(new Point[n]);
    for (int i = 0; i < n; ++i) {
        if (random) {
            points[i].reach = rng.uniformInt(static_cast<unsigned long>(n >> 1));
        } else {
            points[i].reach = i;
        }
    }
    return points;
}

} // namespace

// Tests add() and pop() methods of SeedList class
BOOST_AUTO_TEST_CASE(AddPopBasic) {
    int n = 128;
    // construct points with strictly increasing reachability distance
    boost::shared_array<Point> points = makePoints(n, false);
    SeedList sl(points.get(), n);
    BOOST_CHECK(sl.empty());
    BOOST_CHECK_EQUAL(sl.capacity(), n);
    BOOST_CHECK_EQUAL(sl.pop(), -1);
    BOOST_CHECK(sl.checkInvariants());
    sl.add(0);
    BOOST_CHECK(sl.checkInvariants());
    BOOST_CHECK_EQUAL(sl.pop(), 0);
    BOOST_CHECK_EQUAL(sl.size(), 0);
    BOOST_CHECK(sl.checkInvariants());
    sl.add(n - 1);
    sl.add(0);
    BOOST_CHECK(sl.checkInvariants());
    BOOST_CHECK_EQUAL(sl.pop(), 0);
    BOOST_CHECK(sl.checkInvariants());
    BOOST_CHECK_EQUAL(sl.pop(), n - 1);
    BOOST_CHECK_EQUAL(sl.size(), 0);
    // add points in increasing reachability-distance order
    for (int i = 0; i < n; ++i) {
        sl.add(i);
    }
    BOOST_CHECK_EQUAL(sl.size(), n);
    BOOST_CHECK(sl.checkInvariants());
    // check that points are popped in increasing reachability-distance order
    for (int i = 0; i < n; ++i) {
        BOOST_CHECK_EQUAL(sl.pop(), i);
        BOOST_CHECK(sl.checkInvariants());
    }
    BOOST_CHECK_EQUAL(sl.size(), 0);
    // add points in decreasing reachability-distance order
    for (int i = n - 1; i >= 0; --i) {
        sl.add(i);
    }
    BOOST_CHECK_EQUAL(sl.size(), n);
    BOOST_CHECK(sl.checkInvariants());
    // check that points are popped in increasing reachability-distance order
    for (int i = 0; i < n; ++i) {
        BOOST_CHECK_EQUAL(sl.pop(), i);
        BOOST_CHECK(sl.checkInvariants());
    }
}

// Tests add() and pop() methods of SeedList with randomly ordered inputs
BOOST_AUTO_TEST_CASE(AddPopRandom) {
    int n = 127;
    boost::shared_array<Point> points = makePoints(n, true);
    SeedList sl(points.get(), n);

    for (int i = 0; i < n; ++i) {
        sl.add(i);
        BOOST_CHECK(sl.checkInvariants());
    }
    BOOST_CHECK_EQUAL(sl.size(), n);
    double maxReach = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
        double reach = points[sl.pop()].reach;
        BOOST_CHECK(sl.checkInvariants());
        BOOST_CHECK(reach >= maxReach);
        maxReach = reach;
    }
    BOOST_CHECK_EQUAL(sl.size(), 0);
}

// Tests the update() method of SeedList
BOOST_AUTO_TEST_CASE(Update) {
    int n = 120;
    boost::shared_array<Point> points = makePoints(n, false);
    boost::shared_array<int> order(new int[n]);
    SeedList sl(points.get(), n);
    for (int i = 0; i < n; ++i) {
        sl.add(i);
    }
    for (int i = 0; i < n; ++i) {
        order[i] = sl.pop();
    }
    for (int i = 0; i < n; ++i) {
        sl.add(i);
    }
    for (int i = 0; i < n; ++i) {
        sl.update(i, -points[i].reach);
        BOOST_CHECK(sl.checkInvariants());
    }
    for (int i = 0; i < n; ++i) {
        BOOST_CHECK_EQUAL(sl.pop(), order[n - i - 1]);
    }
}

