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
 

/**
 * @file
 * @brief   Tests for cross-match (distance based and in-ellipse) algorithms.
 *
 * @ingroup associate
 */

#include <algorithm>
#include <iostream>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MatchTest
#include "boost/test/unit_test.hpp"

#include "lsst/afw/math/Random.h"

#include "lsst/ap/EllipseTypes.h"
#include "lsst/ap/Match.h"
#include "lsst/ap/Object.h"
#include "lsst/ap/Point.h"
#include "lsst/ap/Time.h"
#include "lsst/ap/ZoneTypes.h"

using std::size_t;
using boost::int64_t;
using lsst::afw::math::Random;
using namespace lsst::ap;


namespace {

Random & rng() {
    static Random * generator = 0;
    if (generator == 0) {
        TimeSpec ts;
        ts.systemTime();
        generator = new Random(Random::MT19937, static_cast<unsigned long>(ts.tv_sec + ts.tv_nsec));
        std::clog << "\n"
            << "     /\n"
            << "    | Note: Using random number seed " << generator->getSeed() << "\n"
            << "    |       and algorithm " << generator->getAlgorithmName() << "\n"
            << "     \\\n" << std::endl;
    }
    return *generator;
}

int const MAX_EXPECTED = 9;

// A test point or ellipse that knows what it should match
struct TestDatum {

    int64_t _id;
    Point   _loc;
    double  _smaa;
    double  _smia;
    double  _pa;
    int32_t _expected;
    int32_t _matched;
    int64_t _expectId[MAX_EXPECTED];

    TestDatum() { init(-1, Point()); }

    TestDatum(int64_t const id, Point const & p) { init(id, p); }

    void init(int64_t const id, Point const & p);

    bool matchWasExpected(int64_t const id) const;

    int64_t getId()  const { return _id;       }
    double  getRa()  const { return _loc._ra;  }
    double  getDec() const { return _loc._dec; }
    double  getSemiMinorAxisLength() const { return _smia; }
    double  getSemiMajorAxisLength() const { return _smaa; }
    double  getPositionAngle()       const { return _pa;   }
};


void TestDatum::init(int64_t const id, Point const & p) {
    _id       = id;
    _loc      = p;
    _smaa     = 0.0;
    _smia     = 0.0;
    _pa       = 0.0;
    _expected =  0;
    _matched  =  0;
    for(int i = 0; i < MAX_EXPECTED; ++i) {
        _expectId[i] = -1;
    }
}


bool TestDatum::matchWasExpected(int64_t const id) const {
    for(int i = 0; i < _expected; ++i) {
        if (_expectId[i] == id) {
            return true;
        }
    }
    return false;
}


std::ostream & operator<<(std::ostream & o, TestDatum const & d) {
    o << '(' << d._id << ", " << d._loc._ra << ", " << d._loc._dec;
    if (d._expected > 0) {
        o << ") [expecting to match";
    } else {
        o << ") [expecting no matches";
    }
    for (int i = 0; i < d._expected; ++i) {
        o << ' ' << d._expectId[i];
    }
    o << ']';
    return o;
}


// Verify that each entry of the given vector was matched the expected number of times.
void verifyMatchCount(std::vector<TestDatum> const & data) {
    std::vector<TestDatum>::const_iterator end = data.end();
    for (std::vector<TestDatum>::const_iterator i(data.begin()); i != end; ++i) {
        BOOST_CHECK_MESSAGE(i->_matched == i->_expected, "Got " << i->_matched << " matches for " << *i);
    }
}


// Don't need chunks for testing the match algorithm
struct BogusChunk {
    typedef TestDatum Entry;
};


typedef ZoneEntry<BogusChunk>  Ze;
typedef Ellipse<TestDatum>     Ell;
typedef PassthroughFilter<Ze>  Filt;
typedef PassthroughFilter<Ell> EllFilt;
typedef MatchWithDistance<Ze>  Match;
typedef ZoneIndex<Ze>          Zi;
typedef EllipseList<TestDatum> EllList;

} // end of anonymous namespace


// explicit instantiations
template class ZoneEntry<BogusChunk>;
template class Ellipse<TestDatum>;
template class ZoneEntryArray<Ze>;
template class ZoneIndex<Ze>;
template class EllipseList<TestDatum>;


namespace {

// Scrutinizes incoming match lists and validates them
struct MlProcessor {
    typedef MatchWithDistance<Ze> Match;
    typedef std::vector<Match>::iterator MatchIterator;

    double _matchDist;

    MlProcessor(double const d) : _matchDist(d) {}

    void operator()(Ze & entry, MatchIterator begin, MatchIterator end);
};


void MlProcessor::operator()(Ze & entry, MatchIterator begin, MatchIterator end) {

    TestDatum * p = entry._data;
    BOOST_CHECK_MESSAGE(p->_expected == end - begin, "incorrect number of matches for " << *p <<
        ": got " << (end - begin) << " expected " << p->_expected);
    p->_matched += (end - begin);

    for ( ; begin != end; ++begin) {
        TestDatum * s = (*begin)->_data;
        double d = degrees(begin->_distance);
        BOOST_CHECK_MESSAGE(d <= _matchDist, *s << " matched " << *p <<
            " but is " << d << " (> " << _matchDist << ") degrees away");
        BOOST_CHECK_MESSAGE(p->matchWasExpected(s->_id), "unexpected match " << *s << " for " << *p <<
            " separated by " << s->_loc.distance(p->_loc));
        ++s->_matched;
        BOOST_CHECK_MESSAGE(s->_matched == 1, "expected one match for " << *s << " got " << s->_matched);
        BOOST_CHECK_MESSAGE(s->_expectId[0] == p->_id, "unexpected match " << *s << " for " << *p);
    }
}


// Scrutinizes incoming match pairs and validates them
struct MpProcessor {
    void operator()(Ell & e, Ze & z);
};


void MpProcessor::operator()(Ell & e, Ze & z) {
    TestDatum * p = e._data;
    TestDatum * s = z._data;

    ++p->_matched;
    ++s->_matched;
    BOOST_CHECK_MESSAGE(p->matchWasExpected(s->_id), "unexpected match " << *s << " for " << *p);
    BOOST_CHECK_MESSAGE(s->_matched == 1, "expected one match for " << *s << " got " << s->_matched);
    BOOST_CHECK_MESSAGE(s->_expectId[0] == p->_id, "unexpected match " << *s << " for " << *p);
}


// Scrutinizes incoming ellipse match lists and validates them
struct EmlProcessor {
    typedef Ze * Match;
    typedef std::vector<Match>::iterator MatchIterator;

    void operator()(Ell & entry, MatchIterator begin, MatchIterator end);
};


void EmlProcessor::operator()(Ell & e, MatchIterator begin, MatchIterator end) {

    TestDatum * p = e._data;
    BOOST_CHECK_MESSAGE(p->_expected == end - begin, "incorrect number of matches for " << *p <<
        ": got " << (end - begin) << " expected " << p->_expected);
    p->_matched += (end - begin);

    for ( ; begin != end; ++begin) {
        TestDatum * s = (*begin)->_data;
        BOOST_CHECK_MESSAGE(p->matchWasExpected(s->_id), "unexpected match " << *s << " for " << *p);
        ++s->_matched;
        BOOST_CHECK_MESSAGE(s->_matched == 1, "expected one match for " << *s << " got " << s->_matched);
        BOOST_CHECK_MESSAGE(s->_expectId[0] == p->_id, "unexpected match " << *s << " for " << *p);
    }
}

} // end of anonymous namespace


// explicit instantiations of match routines
namespace lsst { namespace ap {
template size_t distanceMatch<Ze, Ze, Filt, Filt, MlProcessor>(
    Zi &,
    Zi &,
    double const,
    Filt &,
    Filt &,
    MlProcessor &
);

template size_t ellipseMatch<TestDatum, Ze, EllFilt, Filt, MpProcessor>(
    EllList &,
    Zi &,
    EllFilt &,
    Filt &,
    MpProcessor &
);

template size_t ellipseGroupedMatch<TestDatum, Ze, EllFilt, Filt, EmlProcessor>(
    EllList &,
    Zi &,
    EllFilt &,
    Filt &,
    EmlProcessor &
);

}} // namespace lsst::ap


namespace {

// builds test data for distance based matching of points
void buildPoints(
    std::vector<TestDatum> & first,
    std::vector<TestDatum> & second,
    Zi                     & fzi,
    Zi                     & szi,
    double const             decMin,
    double const             decMax,
    double const             radius
) {
    BOOST_REQUIRE(decMin < decMax && decMin >= -90.0 && decMax <= 90.0);
    BOOST_REQUIRE(radius > 0.0 && radius < 0.0667);

    double min = clampDec(decMin - 10.0*radius);
    double max = clampDec(decMax + 10.0*radius);
    fzi.setDecBounds(min, max);
    szi.setDecBounds(min, max);

    int64_t i        = 0;
    double  deltaDec = 4.0*radius;

    // divide the test region into ra/dec boxes
    for (double dec = decMin; dec < decMax; dec += deltaDec) {

        double d1 = std::fabs(dec);
        double d2 = std::fabs(dec + deltaDec);
        double const deltaRa = maxAlpha(4.0*radius, clampDec(d1 < d2 ? d2 : d1));

        for (double ra = 0.0; ra <= 360 - 2.0*deltaRa; ra += deltaRa) {
            // create a random point inside a sub region of each box such that a circle of the
            // given radius centered on that point is guaranteed not to cross the box boundaries -
            // we then generate points we know will match by perturbing the initial point.
            TestDatum d;
            if (ra == 0.0) {
                // make sure ra wrap-around is tested
                d.init(
                    i,
                    Point::random(
                        rng(),
                        360.0 - 0.125*deltaRa,
                                0.125*deltaRa,
                        clampDec(dec + deltaDec*0.38),
                        clampDec(dec + deltaDec*0.62)
                    )
                );
            } else {
                d.init(
                    i,
                    Point::random(
                        rng(),
                        ra + deltaRa*0.38,
                        ra + deltaRa*0.62,
                        clampDec(dec + deltaDec*0.38),
                        clampDec(dec + deltaDec*0.62)
                    )
                );
            }
            ++i;

            int nm = static_cast<int>(rng().uniformInt(MAX_EXPECTED + 1));
            assert(nm <= MAX_EXPECTED);
            // generate between 0 and 9 matches for the test point just created
            for (int j = 0; j < nm; ++j) {
                TestDatum m(i, d._loc);
                m._loc.perturb(rng(), radius);
                ++i;
                double const dist = d._loc.distance(m._loc);
                if (dist < 0.95*radius) {
                    d._expectId[d._expected] = m._id;
                    ++d._expected;
                    m._expected = 1;
                    m._expectId[0] = d._id;
                    second.push_back(m);
                }
                else if (dist > 1.05*radius && dist < 1.45*radius) {
                    second.push_back(m);
                }
            }
            first.push_back(d);
        }
    }

    for (size_t i = 0; i < first.size(); ++i) {
        fzi.insert(first[i].getRa(), first[i].getDec(), &first[i], 0, 0);
    }
    for (size_t i = 0; i < second.size(); ++i) {
        szi.insert(second[i].getRa(), second[i].getDec(), &second[i], 0, 0);
    }

    // sort the indexes
    fzi.sort();
    szi.sort();
}


// builds test data for in-ellipse matching
void buildEllipsesAndPoints(
    std::vector<TestDatum> & first,
    std::vector<TestDatum> & second,
    EllList                & ells,
    Zi                     & szi,
    double const             decMin,
    double const             decMax,
    double const             radius
) {
    BOOST_REQUIRE(decMin < decMax && decMin >= -90.0 && decMax <= 90.0);
    BOOST_REQUIRE(radius > 0.0 && radius < 0.0667);

    double min = clampDec(decMin - 10.0*radius);
    double max = clampDec(decMax + 10.0*radius);
    szi.setDecBounds(min, max);

    int64_t i        = 0;
    double  deltaDec = 4.0*radius;

    // divide the test region into ra/dec boxes
    for (double dec = decMin; dec < decMax; dec += deltaDec) {

        double d1 = std::fabs(dec);
        double d2 = std::fabs(dec + deltaDec);
        double const deltaRa = maxAlpha(4.0*radius, clampDec(d1 < d2 ? d2 : d1));

        for (double ra = 0.0; ra <= 360 - 2.0*deltaRa; ra += deltaRa) {
            // create a random point inside a sub region of each box such that a circle of the
            // given radius centered on that point is guaranteed not to cross the box boundaries -
            // we then generate points we know will match by perturbing the initial point.
            TestDatum d;
            if (ra == 0.0) {
                // make sure ra wrap-around is tested
                d.init(
                    i,
                    Point::random(
                        rng(),
                        360.0 - 0.125*deltaRa,
                                0.125*deltaRa,
                        clampDec(dec + deltaDec*0.38),
                        clampDec(dec + deltaDec*0.62)
                    )
                );
            } else {
                d.init(
                    i,
                    Point::random(
                        rng(),
                        ra + deltaRa*0.38,
                        ra + deltaRa*0.62,
                        clampDec(dec + deltaDec*0.38),
                        clampDec(dec + deltaDec*0.62)
                    )
                );
            }
            ++i;
            double tmp = 0.5*std::fabs(rng().gaussian()) + 0.5;
            d._smaa = (tmp >= 1.0 ? 1.0 : tmp)*radius;
            d._smia = rng().flat(0.1, 0.85)*d._smaa;
            d._pa   = rng().flat(0.0, 180.0);

            int nm = static_cast<int>(rng().uniformInt(MAX_EXPECTED + 1));
            assert(nm <= MAX_EXPECTED);
            // generate between 0 and 9 matches for the test ellipses just created
            for (int j = 0; j < nm; ++j) {
                TestDatum m(i, d._loc);
                bool north = rng().uniform() <= 0.5;
                // perturb some amount along major/minor axis, allowing to easily check
                // whether or not the perturbed point should match
                m._loc.perturb(rng(), radius, d._pa + (north ? 0.0 : 90.0));
                ++i;
                double const dist = d._loc.distance(m._loc);
                double cmp = north ? d._smaa : d._smia;
                if (dist < 0.9*cmp) {
                    d._expectId[d._expected] = m._id;
                    ++d._expected;
                    m._expected = 1;
                    m._expectId[0] = d._id;
                    second.push_back(m);
                }
                else if (dist > 1.1*cmp && dist < 1.45*cmp) {
                    second.push_back(m);
                }
            }
            // convert smia/smaa to arc-seconds
            d._smaa *= 3600.0;
            d._smia *= 3600.0;
            first.push_back(d);
        }
    }

    ells.reserve(first.size());
    for (size_t i = 0; i < first.size(); ++i) {
        ells.push_back(Ell(first[i]));
    }
    for (size_t i = 0; i < second.size(); ++i) {
        szi.insert(second[i].getRa(), second[i].getDec(), &second[i], 0, 0);
    }
    szi.sort();
}

} // end of anonymous namespace


BOOST_AUTO_TEST_CASE(distanceMatchTest1) {
    BOOST_TEST_MESSAGE("    - Distance match test, near north pole");
    Zi fzi(60, 60, 1024);
    Zi szi(60, 60, 1024);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    first.reserve(65536);
    second.reserve(65536);
    double const rad = 0.01666666667;
    buildPoints(first, second, fzi, szi, 88.0, 90.0, rad);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " points to " << second.size() << " points ...");
    MlProcessor mlp(rad);
    Filt        f;
    Stopwatch   watch(true);
    size_t      nm = distanceMatch(fzi, szi, rad, f, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(distanceMatchTest2) {
    BOOST_TEST_MESSAGE("    - Distance match test, near south pole");
    Zi fzi(120, 60, 1024);
    Zi szi(30, 45, 1024);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    first.reserve(65536);
    second.reserve(65536);
    double const rad = 0.01666666667;
    buildPoints(first, second, fzi, szi, -90.0, -88.0, rad);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " points to " << second.size() << " points ...");
    MlProcessor mlp(rad);
    Filt        f;
    Stopwatch   watch(true);
    size_t      nm = distanceMatch(fzi, szi, rad, f, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(distanceMatchTest3) {
    BOOST_TEST_MESSAGE("    - Distance match test, near equator");
    Zi fzi(59, 60, 1024);
    Zi szi(44, 60, 1024);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    double const rad = 0.05;
    first.reserve(65536);
    second.reserve(65536);
    buildPoints(first, second, fzi, szi, -1.0, 1.0, rad);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " points to " << second.size() << " points ...");
    MlProcessor mlp(rad);
    Filt        f;
    Stopwatch   watch(true);
    size_t      nm = distanceMatch(fzi, szi, rad, f, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseMatchTest1) {
    BOOST_TEST_MESSAGE("    - Ellipse match test, near north pole");
    Zi                     szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllList                ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, 89.0, 90.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    MpProcessor mpp;
    EllFilt     ef;
    Filt        f;
    Stopwatch   watch(true);
    size_t      nm = ellipseMatch(ells, szi, ef, f, mpp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseMatchTest2) {
    BOOST_TEST_MESSAGE("    - Ellipse match test, near south pole");
    Zi                     szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllList                ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -90.0, -89.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    MpProcessor mpp;
    EllFilt     ef;
    Filt        f;
    Stopwatch   watch(true);
    size_t      nm = ellipseMatch(ells, szi, ef, f, mpp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseMatchTest3) {
    BOOST_TEST_MESSAGE("    - Ellipse match test, near equator");
    Zi                     szi(240, 60, 512);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllList                ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -0.5, 0.5, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    MpProcessor mpp;
    EllFilt     ef;
    Filt        f;
    Stopwatch   watch(true);
    size_t      nm = ellipseMatch(ells, szi, ef, f, mpp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseGroupedMatchTest1) {
    BOOST_TEST_MESSAGE("    - Ellipse grouped match test, near north pole");
    Zi                     szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllList                ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, 89.0, 90.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    EmlProcessor mlp;
    EllFilt      ef;
    Filt         f;
    Stopwatch    watch(true);
    size_t       nm = ellipseGroupedMatch(ells, szi, ef, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseGroupedMatchTest2) {
    BOOST_TEST_MESSAGE("    - Ellipse grouped match test, near south pole");
    Zi                     szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllList                ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -90.0, -89.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    EmlProcessor mlp;
    EllFilt      ef;
    Filt         f;
    Stopwatch    watch(true);
    size_t       nm = ellipseGroupedMatch(ells, szi, ef, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseGroupedMatchTest3) {
    BOOST_TEST_MESSAGE("    - Ellipse grouped match test, near equator");
    Zi                     szi(240, 60, 512);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllList                ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -0.5, 0.5, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    EmlProcessor mlp;
    EllFilt      ef;
    Filt         f;
    Stopwatch    watch(true);
    size_t       nm = ellipseGroupedMatch(ells, szi, ef, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}

