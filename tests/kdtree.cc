#include <algorithm>
#include <iostream>

#include "boost/shared_array.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE KDTree
#include "boost/test/unit_test.hpp"

#include "lsst/afw/math/Random.h"
#include "lsst/ap/Common.h"
#include "lsst/ap/Point.h"
#include "lsst/ap/SpatialUtil.h"
#include "lsst/ap/Time.h"
#include "lsst/ap/cluster/optics/KDTree.cc"
#include "lsst/ap/cluster/optics/Metrics.h"


namespace ap = lsst::ap;
namespace optics = lsst::ap::cluster::optics;

using lsst::afw::math::Random;


namespace {

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

typedef optics::Point<1, void> Point1;
typedef optics::Point<3, MatchOracle> Point;
typedef optics::KDTree<3, MatchOracle> KDTree;

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
// TODO: change to use radians and afw Point/Coord classes when
// rewriting nightly AP code.
void makePoints(std::vector<Point> & points,
                std::vector<MatchOracle> & queryPoints,
                double radiusDeg)
{
    static Random rng(Random::MT19937);
    double minDec = -90.0;
    double maxDec = 90.0;
    double deltaDec = 4.0 * radiusDeg;
    int i = 0;
    // divide the unit sphere into longitude/latitude angle boxes
    for (double dec = minDec; dec < maxDec; dec += deltaDec) {
        double d1 = std::fabs(dec);
        double d2 = std::fabs(dec + deltaDec);
        double const deltaRa = ap::maxAlpha(
            4.0 * radiusDeg, ap::clampDec(std::max(d1, d2)));
        for (double ra = 0.0; ra < 360.0 - deltaRa; ra += deltaRa) {
            // create a random query point inside a sub region of each box
            // such that a circle of the given radius centered on that point
            // is guaranteed not to cross the box boundaries
            MatchOracle m;
            ap::Point qp = ap::Point::random(rng,
                                             ra + deltaRa * 0.38,
                                             ra + deltaRa * 0.62,
                                             ap::clampDec(dec + deltaDec * 0.38),
                                             ap::clampDec(dec + deltaDec * 0.62));
            double theta = ap::RADIANS_PER_DEGREE * qp._ra;
            double phi = ap::RADIANS_PER_DEGREE * qp._dec;
            m.v[0] = std::cos(theta) * std::cos(phi);
            m.v[1] = std::sin(theta) * std::cos(phi);
            m.v[2] = std::sin(phi);
            int nm = static_cast<int>(rng.uniformInt(MatchOracle::MAX_EXPECTED));

            // generate matches for the MatchOracle
            for (int j = 0; j < nm; ++j) {
                ap::Point kdp(qp);
                kdp.perturb(rng, radiusDeg);
                double const dist = kdp.distance(qp);
                if (dist < 1.45 * radiusDeg) {
                    Point p;
                    theta = ap::RADIANS_PER_DEGREE * kdp._ra;
                    phi = ap::RADIANS_PER_DEGREE * kdp._dec;
                    p.coords[0] = std::cos(theta) * std::cos(phi);
                    p.coords[1] = std::sin(theta) * std::cos(phi);
                    p.coords[2] = std::sin(phi);
                    if (dist < 0.999999999 * radiusDeg) {
                        // p is a match - store it and remember its insertion 
                        // index (the kdtree will reorder points)
                        m.expected[m.numExpected++] = i;
                        p.state = i++;
                        points.push_back(p);
                    } else if (dist > 1.0000000001 * radiusDeg) {
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

// tests the median finding algorithm
BOOST_AUTO_TEST_CASE(Median) {
    static Random rng(Random::MT19937);
    // arrays containing 1 coordinate value
    for (int n = 1; n < 50; ++n) {
        boost::shared_array<Point1> points(new Point1[n]);
        for (int i = 0; i < n; ++i) {
            points[i].coords[0] = 1.0;
        }
        BOOST_CHECK_EQUAL(optics::median(points.get(), n, 0), n >> 1);
    }
    // arrays containing distinct coordinate values
    for (int n = 1; n < 100; ++n) {
        boost::shared_array<Point1> points(new Point1[n]);
        for (int i = 0; i < n; ++i) {
            points[i].coords[0] = i;
        }
        BOOST_CHECK_EQUAL(optics::median(points.get(), n, 0), n >> 1);
        double medianValue = points[n >> 1].coords[0];
        std::reverse(points.get(), points.get() + n);
        int m = optics::median(points.get(), n, 0);
        BOOST_CHECK_EQUAL(points[m].coords[0], medianValue);
        shuffle(points.get(), points.get() + n);
        m = optics::median(points.get(), n, 0);
        BOOST_CHECK_EQUAL(points[m].coords[0], medianValue);
    }
    // arrays containing repeated values
    for (int n = 2; n < 50; ++n) {
        std::vector<Point1> points;
        for (int i = 0; i < n; ++i) {
            unsigned long nr = rng.uniformInt(10u) + 1;
            for (unsigned long j = 0; j < nr; ++j) {
                Point1 p;
                p.coords[0] = i;
                points.push_back(p);
            }
        }
        int sz = static_cast<int>(points.size());
        int medianIndex = sz >> 1;
        double medianValue = points[medianIndex].coords[0];
        BOOST_CHECK_EQUAL(optics::median(&points[0], sz, 0), medianIndex);
        std::reverse(points.begin(), points.end());
        int m = optics::median(&points[0], sz, 0);
        BOOST_CHECK_EQUAL(points[m].coords[0], medianValue);
        shuffle(points.begin(), points.end());
        m = optics::median(&points[0], sz, 0);
        BOOST_CHECK_EQUAL(points[m].coords[0], medianValue);
    }
}

// tests kdtree range search
BOOST_AUTO_TEST_CASE(RangeQuery) {
    typedef std::vector<MatchOracle>::const_iterator Iter;

    double const radiusDeg = 0.5;
    double const radiusRad = ap::RADIANS_PER_DEGREE * radiusDeg;
    double d = 2.0 * std::sin(0.5 * radiusRad);
    d = d * d;

    optics::SquaredEuclidianDistanceOverSphere metric;
    std::vector<Point> points;
    std::vector<MatchOracle> queryPoints;
    makePoints(points, queryPoints, radiusDeg);
    BOOST_TEST_MESSAGE("Number of query points: " << queryPoints.size());
    BOOST_TEST_MESSAGE("Number of points: " << points.size());
    ap::Stopwatch watch(true); 
    KDTree tree(&points[0], static_cast<int>(points.size()), 32, 0.0);
    watch.stop();
    BOOST_TEST_MESSAGE("Built tree in " << watch);
    watch.start();
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
    watch.stop();
    BOOST_TEST_MESSAGE("Validated range query for each query point in " << watch);
}

