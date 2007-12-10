// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   MatchTest.cc
//! \brief  Tests for cross-match (distance based and in-ellipse) algorithms.
//
//##====----------------                                ----------------====##/

#include <algorithm>
#include <iostream>
#include <vector>

#include <boost/version.hpp>
#if BOOST_VERSION < 103400
#   include <boost/test/auto_unit_test.hpp>
#   define BOOST_TEST_MESSAGE BOOST_MESSAGE
#else
#   include <boost/test/unit_test.hpp>
#endif

#include <lsst/ap/EllipseTypes.h>
#include <lsst/ap/Match.h>
#include <lsst/ap/Object.h>
#include <lsst/ap/Random.h>
#include <lsst/ap/Time.h>
#include <lsst/ap/ZoneTypes.h>

using namespace lsst::ap;


namespace {

int32_t const MAX_EXPECTED = 9;

// A test point/ellipse that knows what it should match
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
    _id          = id;
    _loc         = p;
    _smaa        = 0.0;
    _smia        = 0.0;
    _pa          = 0.0;
    _expected    =  0;
    _matched     =  0;
    for(int32_t i = 0; i < MAX_EXPECTED; ++i) {
        _expectId[i] = -1;
    }
}


bool TestDatum::matchWasExpected(int64_t const id) const {
    for(int32_t i = 0; i < _expected; ++i) {
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
    for (int32_t i = 0; i < d._expected; ++i) {
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
    typedef TestDatum EntryType;
};


// explicit instantiations and helpful typedefs
template class ZoneEntry<BogusChunk>;
template class Ellipse<TestDatum>;

typedef ZoneEntry<BogusChunk>      ZeType;
typedef Ellipse<TestDatum>         EllType;
typedef PassthroughFilter<ZeType>  FiltType;
typedef PassthroughFilter<EllType> EllFiltType;
typedef MatchWithDistance<ZeType>  MatchType;

template class Zone<ZeType>;
template class ZoneIndex<ZeType>;
template class EllipseList<TestDatum>;

typedef ZoneIndex<ZeType>      ZiType;
typedef EllipseList<TestDatum> EllListType;


// Scrutinizes incoming match lists and validates them
struct MlProcessor {
    typedef MatchWithDistance<ZeType> MatchType;
    typedef std::vector<MatchType>::iterator MatchIteratorType;

    double _matchDist;

    MlProcessor(double const d) : _matchDist(d) {}

    void operator()(ZeType & entry, MatchIteratorType begin, MatchIteratorType end);
};


void MlProcessor::operator()(ZeType & entry, MatchIteratorType begin, MatchIteratorType end) {

    TestDatum * p = entry._data;
    BOOST_CHECK_MESSAGE(p->_expected == end - begin, "incorrect number of matches for " << *p <<
        ": got " << (end - begin) << " expected " << p->_expected);
    p->_matched += (end - begin);

    for ( ; begin != end; ++begin) {
        TestDatum * s = (*begin)->_data;
        BOOST_CHECK_MESSAGE((*begin)->_distance <= _matchDist, *s << " matched " << *p <<
            " but is " << (*begin)->_distance << " (> " << _matchDist << ") degrees away");
        BOOST_CHECK_MESSAGE(p->matchWasExpected(s->_id), "unexpected match " << *s << " for " << *p <<
            " separated by " << s->_loc.distance(p->_loc));
        ++s->_matched;
        BOOST_CHECK_MESSAGE(s->_matched == 1, "expected one match for " << *s << " got " << s->_matched);
        BOOST_CHECK_MESSAGE(s->_expectId[0] == p->_id, "unexpected match " << *s << " for " << *p);
    }
}


// Scrutinizes incoming match pairs and validates them
struct MpProcessor {
    void operator()(EllType & e, ZeType & z);
};


void MpProcessor::operator()(EllType & e, ZeType & z) {
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
    typedef ZeType * MatchType;
    typedef std::vector<MatchType>::iterator MatchIteratorType;

    void operator()(EllType & entry, MatchIteratorType begin, MatchIteratorType end);
};


void EmlProcessor::operator()(EllType & e, MatchIteratorType begin, MatchIteratorType end) {

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


// explicit instantiations of match routines

template size_t distanceMatch<ZeType, ZeType, FiltType, FiltType, MlProcessor>(
    ZiType &,
    ZiType &,
    double const,
    FiltType &,
    FiltType &,
    MlProcessor &
);

template size_t ellipseMatch<TestDatum, ZeType, EllFiltType, FiltType, MpProcessor>(
    EllListType &,
    ZiType &,
    EllFiltType &,
    FiltType &,
    MpProcessor &
);

template size_t ellipseGroupedMatch<TestDatum, ZeType, EllFiltType, FiltType, EmlProcessor>(
    EllListType &,
    ZiType &,
    EllFiltType &,
    FiltType &,
    EmlProcessor &
);


// builds test data for distance based matching of points
void buildPoints(
    std::vector<TestDatum> & first,
    std::vector<TestDatum> & second,
    ZiType &                 fzi,
    ZiType &                 szi,
    double const             decMin,
    double const             decMax,
    double const             radius
) {
    BOOST_REQUIRE(decMin < decMax && decMin >= -90.0 && decMax <= 90.0);
    BOOST_REQUIRE(radius > 0.0 && radius < 0.0667);

    initRandom();

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
                        ra + deltaRa*0.38,
                        ra + deltaRa*0.62,
                        clampDec(dec + deltaDec*0.38),
                        clampDec(dec + deltaDec*0.62)
                    )
                );
            }
            ++i;

            int32_t nm = static_cast<int32_t>(std::floor(uniformRandom()*(MAX_EXPECTED + 1)));
            assert(nm <= MAX_EXPECTED);
            // generate between 0 and 9 matches for the test point just created
            for (int32_t j = 0; j < nm; ++j) {
                TestDatum m(i, d._loc);
                m._loc.perturb(radius);
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
        fzi.insert(&first[i], 0, 0);
    }
    for (size_t i = 0; i < second.size(); ++i) {
        szi.insert(&second[i], 0, 0);
    }

    // sort the indexes
    fzi.sort();
    szi.sort();
}


// builds test data for in-ellipse matching
void buildEllipsesAndPoints(
    std::vector<TestDatum> & first,
    std::vector<TestDatum> & second,
    EllListType            & ells,
    ZiType &                 szi,
    double const             decMin,
    double const             decMax,
    double const             radius
) {
    BOOST_REQUIRE(decMin < decMax && decMin >= -90.0 && decMax <= 90.0);
    BOOST_REQUIRE(radius > 0.0 && radius < 0.0667);

    initRandom();

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
                        ra + deltaRa*0.38,
                        ra + deltaRa*0.62,
                        clampDec(dec + deltaDec*0.38),
                        clampDec(dec + deltaDec*0.62)
                    )
                );
            }
            ++i;
            double tmp = 0.5*std::fabs(normalRandom()) + 0.5;
            d._smaa = (tmp >= 1.0 ? 1.0 : tmp)*radius;
            d._smia = (0.75*uniformRandom() + 0.1)*d._smaa;
            d._pa   = uniformRandom()*180.0;

            int32_t nm = static_cast<int32_t>(std::floor(uniformRandom()*(MAX_EXPECTED + 1)));
            assert(nm <= MAX_EXPECTED);
            // generate between 0 and 9 matches for the test ellipses just created
            for (int32_t j = 0; j < nm; ++j) {
                TestDatum m(i, d._loc);
                bool north = coinToss(0.5);
                // perturb some amount along major/minor axis, allowing to easily check
                // whether or not the perturbed point should match
                m._loc.perturb(radius, d._pa + (north ? 0.0 : 90.0));
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
        ells.push_back(first[i]);
    }
    for (size_t i = 0; i < second.size(); ++i) {
        szi.insert(&second[i], 0, 0);
    }
    szi.sort();
}

} // end of anonymous namespace


BOOST_AUTO_TEST_CASE(distanceMatchTest1) {
    BOOST_TEST_MESSAGE("    - Distance match test, near north pole");
    ZiType fzi(60, 60, 1024);
    ZiType szi(60, 60, 1024);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    first.reserve(65536);
    second.reserve(65536);
    double const rad = 0.01666666667;
    buildPoints(first, second, fzi, szi, 88.0, 90.0, rad);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " points to " << second.size() << " points ...");
    MlProcessor mlp(rad);
    FiltType    f;
    Stopwatch   watch(true);
    size_t      nm = distanceMatch(fzi, szi, rad, f, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(distanceMatchTest2) {
    BOOST_TEST_MESSAGE("    - Distance match test, near south pole");
    ZiType fzi(120, 60, 1024);
    ZiType szi(30, 45, 1024);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    first.reserve(65536);
    second.reserve(65536);
    double const rad = 0.01666666667;
    buildPoints(first, second, fzi, szi, -90.0, -88.0, rad);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " points to " << second.size() << " points ...");
    MlProcessor mlp(rad);
    FiltType    f;
    Stopwatch   watch(true);
    size_t      nm = distanceMatch(fzi, szi, rad, f, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(distanceMatchTest3) {
    BOOST_TEST_MESSAGE("    - Distance match test, near equator");
    ZiType fzi(59, 60, 1024);
    ZiType szi(44, 60, 1024);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    double const rad = 0.05;
    first.reserve(65536);
    second.reserve(65536);
    buildPoints(first, second, fzi, szi, -1.0, 1.0, rad);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " points to " << second.size() << " points ...");
    MlProcessor mlp(rad);
    FiltType    f;
    Stopwatch   watch(true);
    size_t      nm = distanceMatch(fzi, szi, rad, f, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseMatchTest1) {
    BOOST_TEST_MESSAGE("    - Ellipse match test, near north pole");
    ZiType                 szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllListType            ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, 89.0, 90.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    MpProcessor mpp;
    EllFiltType ef;
    FiltType    f;
    Stopwatch   watch(true);
    size_t      nm = ellipseMatch(ells, szi, ef, f, mpp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseMatchTest2) {
    BOOST_TEST_MESSAGE("    - Ellipse match test, near south pole");
    ZiType                 szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllListType            ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -90.0, -89.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    MpProcessor mpp;
    EllFiltType ef;
    FiltType    f;
    Stopwatch   watch(true);
    size_t      nm = ellipseMatch(ells, szi, ef, f, mpp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseMatchTest3) {
    BOOST_TEST_MESSAGE("    - Ellipse match test, near equator");
    ZiType                 szi(240, 60, 512);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllListType            ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -0.5, 0.5, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    MpProcessor mpp;
    EllFiltType ef;
    FiltType    f;
    Stopwatch   watch(true);
    size_t      nm = ellipseMatch(ells, szi, ef, f, mpp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseGroupedMatchTest1) {
    BOOST_TEST_MESSAGE("    - Ellipse grouped match test, near north pole");
    ZiType                 szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllListType            ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, 89.0, 90.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    EmlProcessor mlp;
    EllFiltType  ef;
    FiltType     f;
    Stopwatch    watch(true);
    size_t       nm = ellipseGroupedMatch(ells, szi, ef, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseGroupedMatchTest2) {
    BOOST_TEST_MESSAGE("    - Ellipse grouped match test, near south pole");
    ZiType                 szi(240, 60, 128);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllListType            ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -90.0, -89.0, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    EmlProcessor mlp;
    EllFiltType  ef;
    FiltType     f;
    Stopwatch    watch(true);
    size_t       nm = ellipseGroupedMatch(ells, szi, ef, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}


BOOST_AUTO_TEST_CASE(ellipseGroupedMatchTest3) {
    BOOST_TEST_MESSAGE("    - Ellipse grouped match test, near equator");
    ZiType                 szi(240, 60, 512);
    std::vector<TestDatum> first;
    std::vector<TestDatum> second;
    EllListType            ells;
    first.reserve(65536);
    second.reserve(65536);
    buildEllipsesAndPoints(first, second, ells, szi, -0.5, 0.5, 0.01666666667);
    BOOST_TEST_MESSAGE("      matching " << first.size() << " ellipses to " << second.size() << " points ...");
    EmlProcessor mlp;
    EllFiltType  ef;
    FiltType     f;
    Stopwatch    watch(true);
    size_t       nm = ellipseGroupedMatch(ells, szi, ef, f, mlp);
    watch.stop();
    BOOST_TEST_MESSAGE("      found and validated " << nm << " match pairs in " << watch);
    verifyMatchCount(first);
    verifyMatchCount(second);
}

