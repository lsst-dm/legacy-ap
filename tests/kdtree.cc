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
 
#include <algorithm>
#include <iostream>

#include "boost/shared_array.hpp"
#include "boost/timer.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE KDTree
#include "boost/test/unit_test.hpp"

#include "lsst/afw/math/Random.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/ap/utils/SpatialUtils.h"
#include "lsst/ap/cluster/optics/KDTree.cc"
#include "lsst/ap/cluster/optics/Metrics.h"


namespace utils = lsst::ap::utils;
namespace optics = lsst::ap::cluster::optics;
namespace afwGeom = lsst::afw::geom;

using lsst::afw::math::Random;
using lsst::afw::geom::Angle;
using lsst::afw::geom::degrees;
using lsst::afw::geom::radians;
using lsst::afw::geom::HALFPI;
using lsst::afw::geom::TWOPI;
using lsst::afw::coord::IcrsCoord;


namespace {

// Pick a declination uniformly at random in the given range
Angle randomDec(Random & rng, Angle decMin, Angle decMax) {
    double z = rng.flat(std::sin(decMin.asRadians()), std::sin(decMax.asRadians()));
    Angle res = std::asin(z) * radians;
    if (res < decMin) {
        return decMin;
    } else if (res > decMax) {
        return decMax;
    }
    return res;
}

// Pick a point uniformly at random in the specified box.
IcrsCoord const random(Random & rng, Angle raMin, Angle raMax, Angle decMin, Angle decMax) {
    assert(raMin < raMax);
    Angle ra = rng.flat(raMin.asRadians(), raMax.asRadians()) * radians;
    return IcrsCoord(ra, randomDec(rng, decMin, decMax));
}

// Perturb coords according to a normal distribution in the direction given by the specified position angle.
IcrsCoord const perturb(Random & rng, IcrsCoord const & coords, Angle sigma, Angle pa) {
    double sra = std::sin(coords.getLongitude().asRadians());
    double cra = std::cos(coords.getLongitude().asRadians());
    double sde = std::sin(coords.getLatitude().asRadians());
    double cde = std::cos(coords.getLatitude().asRadians());
    double spa = std::sin(pa.asRadians());
    double cpa = std::cos(pa.asRadians());

    // 3-vector corresponding to coords
    Eigen::Vector3d p(cra*cde, sra*cde, sde);
    // north vector tangential to p
    Eigen::Vector3d n(-cra*sde, -sra*sde, cde);
    // east vector tangential to p
    Eigen::Vector3d e(-sra, cra, 0.0);
    // rotate north vector at V by minus position angle
    Eigen::Vector3d nt = spa*e + cpa*n;
    // get magnitude of gaussian perturbation
    double m = rng.gaussian() * sigma.asRadians();
    // return the perturbed position
    return utils::cartesianToIcrs(p*std::cos(m)+ nt*std::sin(m));
}

// Perturb coords according to a normal distribution.
IcrsCoord const perturb(Random & rng, IcrsCoord const & coords, Angle sigma) {
    return perturb(rng, coords, sigma, rng.uniform() * TWOPI * radians);
}


struct MatchOracle {
    static int const MAX_EXPECTED = 57;

    Eigen::Vector3d v;
    int numExpected;
    int expected[MAX_EXPECTED];

    MatchOracle() : numExpected(0) {
        for (int j = 0; j < MAX_EXPECTED; ++j) {
            expected[j] = -1;
        }
    }

    bool isExpected(int i) const {
        return std::binary_search(expected, expected + numExpected, i);
    }
};

typedef optics::Point<3, MatchOracle *> Point;
typedef optics::KDTree<3, MatchOracle *> KDTree;

// randomly shuffles an input sequence
template <typename RandomAccessIterT>
void shuffle(RandomAccessIterT begin, RandomAccessIterT end) {
    static Random rng(Random::MT19937);
    using std::swap;
    unsigned long n = static_cast<unsigned long>(end - begin);

    for (unsigned long i = 0; i < n; ++i) {
        RandomAccessIterT j = begin + rng.uniformInt(n);
        RandomAccessIterT k = begin + rng.uniformInt(n);
        swap(*j, *k);
    }
}

// creates test points for KDTree range search
void makePoints(std::vector<Point> & points,
                std::vector<MatchOracle> & queryPoints,
                Angle radius)
{
    static Random rng(Random::MT19937);
    Angle minDec = -HALFPI * radians;
    Angle maxDec =  HALFPI * radians;
    Angle deltaDec = 4.0 * radius;
    int i = 0;
    // divide the unit sphere into longitude/latitude angle boxes
    for (Angle dec = minDec; dec < maxDec; dec += deltaDec) {
        Angle d1 = std::fabs(dec) * radians;
        Angle d2 = std::fabs(dec + deltaDec) * radians;
        Angle const deltaRa = utils::maxAlpha(
            4.0 * radius, utils::clampPhi(std::max(d1, d2)));
        for (Angle ra = 0.0 * radians; ra < TWOPI * radians - deltaRa; ra += deltaRa) {
            // create a random query point inside a sub region of each box
            // such that a circle of the given radius centered on that point
            // is guaranteed not to cross the box boundaries
            MatchOracle m;
            IcrsCoord qp = random(
                rng,
                ra + deltaRa * 0.38,
                ra + deltaRa * 0.62,
                utils::clampPhi(dec + deltaDec * 0.38),
                utils::clampPhi(dec + deltaDec * 0.62)
            );
            m.v = qp.getVector().asEigen();
            int nm = static_cast<int>(rng.uniformInt(MatchOracle::MAX_EXPECTED));

            // generate matches for the MatchOracle
            for (int j = 0; j < nm; ++j) {
                IcrsCoord pp = perturb(rng, qp, radius);
                Angle const dist = pp.angularSeparation(qp);
                if (dist < 1.45 * radius) {
                    Point p;
                    p.coords = pp.getVector().asEigen();
                    if (dist < 0.999999999 * radius) {
                        // p is a match - store it and remember its insertion 
                        // index (the kdtree will reorder points)
                        m.expected[m.numExpected++] = i;
                        p.state = i++;
                        points.push_back(p);
                    } else if (dist > 1.0000000001 * radius) {
                        // p does not match - store it and save its insertion
                        // index
                        p.state = i++;
                        points.push_back(p);
                    }
                }
            }
            queryPoints.push_back(m);
        }
    }
}

} // namespace

// tests kdtree range search
BOOST_AUTO_TEST_CASE(RangeQuery) {
    typedef std::vector<MatchOracle>::const_iterator Iter;

    Angle const radius = 0.5 * degrees;
    double d = 2.0 * std::sin(0.5 * radius.asRadians());
    d = d * d;

    optics::SquaredEuclidianDistanceOverSphere metric;
    std::vector<Point> points;
    std::vector<MatchOracle> queryPoints;
    makePoints(points, queryPoints, radius);
    BOOST_TEST_MESSAGE("Number of query points: " << queryPoints.size());
    BOOST_TEST_MESSAGE("Number of points: " << points.size());
    boost::timer watch; 
    KDTree tree(&points[0], static_cast<int>(points.size()), 32, 0.0);
    BOOST_TEST_MESSAGE("Built tree in " << watch.elapsed() << " s");
    watch.restart();
    for (Iter i = queryPoints.begin(); i != queryPoints.end(); ++i) {
        MatchOracle const & m = *i;
        int nm = 0;
        int j = tree.inRange(m.v, d, metric);
        while (j != -1) {
            ++nm;
            BOOST_CHECK(m.isExpected(points[j].state));
            BOOST_CHECK(points[j].dist <= d);
            j = points[j].next;
        }
        BOOST_CHECK_EQUAL(m.numExpected, nm);
    }
    BOOST_TEST_MESSAGE("Validated range query for each query point in " <<
                       watch.elapsed() << " s");
}

