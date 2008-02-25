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
#include <map>

#include <boost/version.hpp>
#include <boost/bind.hpp>
#if BOOST_VERSION < 103400
#   include <boost/test/auto_unit_test.hpp>
#   define BOOST_TEST_MESSAGE BOOST_MESSAGE
#else
#   include <boost/test/unit_test.hpp>
#endif

#include <lsst/ap/ChunkManagerImpl.cc>
#include <lsst/ap/Random.h>


using namespace lsst::ap;


struct TestEntry {
    int64_t  _id;
    int      _next;

    TestEntry() : _id(-1), _next(-1) {}

    int64_t getId()          const           { return _id;   }
    int     getNextInChain() const           { return _next; }
    void    setId         (int64_t const id) { _id   = id;   }
    void    setNextInChain(int     const id) { _next = id;   }
};


typedef std::vector<int64_t>               TestIds;
typedef detail::HashedSet<TestEntry, 8192> TestEntrySet;

template class detail::HashedSet<TestEntry, 8192>;


void initTestIds(TestIds & ids, int const n) {
    BOOST_REQUIRE(n > 0);
    initRandom();
    ids.clear();
    ids.reserve(n);
    for (int i = 0; i < n; ++i) {
        ids.push_back(static_cast<int64_t>(std::ldexp(uniformRandom() - 0.5, 62)));
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
    BOOST_CHECK_EQUAL(hs->size(), v.size());
    for (TestIds::const_iterator i = v.begin(); i != end; ++i) {
        BOOST_CHECK_MESSAGE(hs->find(*i) != 0, "lost previously inserted item " << *i);
        int64_t outside = static_cast<int64_t>(std::ldexp(uniformRandom(), 62));
        BOOST_CHECK_MESSAGE(hs->find(outside) == 0, "found item " << outside << " that was never inserted");
    }
    // fill up the set with random entries
    while (hs->space() > 0) {
        hs->insert(static_cast<int64_t>(std::ldexp(uniformRandom() - 0.5, 62)));
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
    BOOST_CHECK_EQUAL(hs->size(), v.size());
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
    BOOST_CHECK_EQUAL(hs->size(), 0u);

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
        BOOST_CHECK_EQUAL(hs->size(), v2.size());
        std::swap(v, v2);
    }

    for (TestIds::const_iterator i = v.begin(); i != v.end(); ++i) {
        BOOST_CHECK_MESSAGE(hs->erase(*i) == true, "failed to erase item " << *i);
    }
    BOOST_CHECK_EQUAL(hs->size(), 0u);
}


BOOST_AUTO_TEST_CASE(hsTest3) {
    BOOST_TEST_MESSAGE("    - HashedSet vs. std::set comparison test");
    // generate test ids
    TestIds v;
    TestIds v2;
    initTestIds(v, 4096);
    initTestIds(v2, 4096);

    // Use STL analogue of HashedSet for reference behaviour
    std::set<int64_t> ref;
    boost::scoped_ptr<TestEntrySet> hs(new TestEntrySet());

    for (TestIds::const_iterator i = v.begin(); i != v.end(); ++i) {
        BOOST_CHECK_MESSAGE(hs->insert(*i) != 0, "ran out of space for new items");
        ref.insert(*i);
    }

    // perform a large number of inserts/erases
    for (int num_passes = 0; num_passes < 50; ++num_passes) {
        for (TestIds::const_iterator i = v2.begin(); i != v2.end(); ++i) {
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
        std::swap(v, v2);
    }

    BOOST_CHECK_EQUAL(hs->size(), ref.size());
    for (std::set<int64_t>::const_iterator i = ref.begin(); i != ref.end(); ++i) {
        BOOST_CHECK_MESSAGE(hs->find(*i) != 0, "std::set and HashedSet disagree");
    }
}

