// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Tests for HashedSet template.
 *
 * @ingroup associate
 */

#include <cmath>
#include <algorithm>
#include <vector>
#include <set>

#include "boost/bind.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HashedSetTest
#include "boost/test/unit_test.hpp"

#include "lsst/afw/math/Random.h"

#include "lsst/ap/ChunkManagerImpl.cc"

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


struct TestEntry {
    int _id;
    int _next;

    TestEntry() : _id(-1), _next(-1) {}

    int getId() const { return _id; }
    int getNextInChain() const { return _next; }
    void setId(int const id) { _id = id; }
    void setNextInChain(int const id) { _next = id; }
};


typedef std::vector<int> TestIds;
typedef detail::HashedSet<TestEntry, 8192> TestEntrySet;

}

template class detail::HashedSet<TestEntry, 8192>;


void initTestIds(TestIds & ids, int const n) {
    BOOST_REQUIRE(n > 0);
    ids.clear();
    ids.reserve(n);
    for (int i = 0; i < n; ++i) {
        int id = static_cast<int>(std::ldexp(rng().flat(-0.5, 0.5), 30));
        ids.push_back(id);
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}


BOOST_AUTO_TEST_CASE(hsTest1) {
    BOOST_TEST_MESSAGE("    - HashedSet find/insert test");
    TestIds v;
    initTestIds(v, 8192);
    boost::scoped_ptr<TestEntrySet> hs(new TestEntrySet());
    TestIds::const_iterator const end = v.end();
    // make sure an empty hash table is properly initialized
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        BOOST_CHECK_MESSAGE(hs->find(*i) == 0, "found item " << *i << " before it was inserted");
    }
    // insert test elements
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        BOOST_CHECK_MESSAGE(hs->find(*i) == 0, "found item " << *i << " before it was inserted");
        BOOST_CHECK_MESSAGE(hs->insert(*i) != 0, "ran out of space for new items");
    }
    // make sure they can be found
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        BOOST_CHECK_MESSAGE(hs->insert(*i) == 0, "Item " << *i << " should already have been inserted");
    }
    // and all have been found
    BOOST_CHECK_EQUAL(hs->size(), static_cast<int>(v.size()));
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        BOOST_CHECK_MESSAGE(hs->find(*i) != 0, "lost previously inserted item " << *i);
        int outside = static_cast<int>(std::ldexp(rng().uniform(), 30));
        BOOST_CHECK_MESSAGE(hs->find(outside) == 0, "found item " << outside << " that was never inserted");
    }
    // fill up the set with random entries
    while (hs->space() > 0) {
        hs->insert(static_cast<int>(std::ldexp(rng().flat(-0.5, 0.5), 30)));
    }
    BOOST_CHECK_MESSAGE(hs->insert(1) == 0, "Inserting into a full HashedSet should fail");
}


BOOST_AUTO_TEST_CASE(hsTest2) {
    BOOST_TEST_MESSAGE("    - HashedSet findOrInsert/erase test");
    // generate test ids and split them across two vectors
    TestIds v;
    initTestIds(v, 8192);
    TestIds::iterator middle = v.begin();
    std::advance(middle, v.size()/2);
    TestIds v2(middle, v.end());
    v.erase(middle, v.end());

    // insert some test elements
    boost::scoped_ptr<TestEntrySet> hs(new TestEntrySet());
    TestIds::const_iterator const end = v.end();
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        std::pair<TestEntry *, bool> p = hs->findOrInsert(*i);
        BOOST_CHECK_MESSAGE(p.first != 0, "ran out of space for new items");
        BOOST_CHECK_MESSAGE(p.second == true, "found pre-existing item with id " << *i);
    }
    BOOST_CHECK_EQUAL(hs->size(), static_cast<int>(v.size()));
    // make sure we can still find them
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        std::pair<TestEntry *, bool> p = hs->findOrInsert(*i);
        BOOST_CHECK_MESSAGE(p.first != 0, "failed to find previously inserted item " << *i);
        BOOST_CHECK_MESSAGE(p.second == false, "failed to find pre-existing item " << *i);
    }
    // erase them all
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        BOOST_CHECK_MESSAGE(hs->erase(*i) == true, "failed to erase item " << *i << ": item not found");
        BOOST_CHECK_MESSAGE(hs->erase(*i) == false, "failed to erase item " << *i);
        BOOST_CHECK_MESSAGE(hs->find(*i) == 0, "found item " << *i << " after it was erased");
    }
    BOOST_CHECK_EQUAL(hs->size(), 0);

    // insert elements of v
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        BOOST_CHECK_MESSAGE(hs->insert(*i) != 0, "ran out of space for new items");
    }
    // perform a a few insert/erase cycles
    for (int num_passes = 0; num_passes < 10; ++num_passes) {
        for (TestIds::const_iterator i = v2.begin(); i != v2.end(); ++i) {
            std::pair<TestEntry *, bool> p = hs->findOrInsert(*i);
            BOOST_CHECK_MESSAGE(p.first != 0, "ran out of space for new items");
            BOOST_CHECK_MESSAGE(p.second == true, "found pre-existing item with id " << *i);
        }
        for (TestIds::const_iterator i = v.begin(); i != v.end(); ++i) {
            BOOST_CHECK_MESSAGE(hs->erase(*i) == true, "failed to erase item " << *i);
        }
        BOOST_CHECK_EQUAL(hs->size(), static_cast<int>(v2.size()));
        std::swap(v, v2);
    }

    for (TestIds::const_iterator i = v.begin(); i != v.end(); ++i) {
        BOOST_CHECK_MESSAGE(hs->erase(*i) == true, "failed to erase item " << *i);
    }
    BOOST_CHECK_EQUAL(hs->size(), 0);
}


BOOST_AUTO_TEST_CASE(hsTest3) {
    BOOST_TEST_MESSAGE("    - HashedSet vs. std::set comparison test");
    // generate test ids
    TestIds v;
    initTestIds(v, 4096);

    // Use STL analogue of HashedSet for reference behaviour
    std::set<int> ref;
    boost::scoped_ptr<TestEntrySet> hs(new TestEntrySet());

    // perform a large number of inserts/erases
    for (int num_passes = 0; num_passes < 50; ++num_passes) {
        for (TestIds::const_iterator i = v.begin(); i != v.end(); ++i) {
            std::pair<TestEntry *, bool> p = hs->findOrInsert(*i);
            BOOST_CHECK_MESSAGE(p.first != 0, "ran out of space for new items");
            BOOST_CHECK_MESSAGE(p.second == (ref.find(*i) == ref.end()),
                "found pre-existing item with id " << *i);
            ref.insert(*i);
        }
        for (TestIds::const_iterator i = v.begin(); i != v.end(); ++i) {
            BOOST_CHECK_MESSAGE(hs->erase(*i) == true, "failed to erase item " << *i);
            ref.erase(*i);
        }
        initTestIds(v, 4096);
    }

    BOOST_CHECK_EQUAL(hs->size(), static_cast<int>(ref.size()));
    for (std::set<int>::const_iterator i = ref.begin(); i != ref.end(); ++i) {
        BOOST_CHECK_MESSAGE(hs->find(*i) != 0, "std::set and HashedSet disagree");
    }
}

