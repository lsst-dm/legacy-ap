// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Tests whether ChunkManager handles successive visits
 *          (including overlapping ones) correctly.
 *
 * @ingroup associate
 */

#include <vector>

#include <boost/version.hpp>
#include <boost/bind.hpp>
#if BOOST_VERSION < 103400
#   include <boost/test/auto_unit_test.hpp>
#   define BOOST_TEST_MESSAGE BOOST_MESSAGE
#else
#   include <boost/test/unit_test.hpp>
#endif

#include <lsst/ap/Common.h>
#include <lsst/ap/Chunk.h>
#include <lsst/ap/ChunkManager.h>
#include <lsst/ap/ScopeGuard.h>
#include <lsst/ap/Time.h>


using namespace lsst::ap;

typedef SharedSimpleObjectChunkManager::SimpleObjectChunk SObjChunk;


BOOST_AUTO_TEST_CASE(disjointVisitsTest) {

    BOOST_TEST_MESSAGE("    - ChunkManager test: sequence of disjoint visits");
    SharedSimpleObjectChunkManager mgr("test");
    // unlink the shared memory object immediately (it remains available until the test process exits)
    SharedSimpleObjectChunkManager::destroyInstance("test");

    // Process a series of non-overlapping visits
    static int64_t const numVisits = 50;
    for (int64_t visitId = 1; visitId < numVisits; ++visitId) {

        ScopeGuard v1(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, visitId - 1, false));
        ScopeGuard v2(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, visitId, false));
        std::vector<SObjChunk> toRead;
        std::vector<SObjChunk> toWaitFor;
        std::vector<int64_t>   chunkIds;

        mgr.registerVisit(visitId);
        chunkIds.push_back(2*visitId);
        chunkIds.push_back(2*visitId + 1);
        mgr.startVisit(toRead, toWaitFor, visitId, chunkIds);
        if (visitId > 1) {
            BOOST_CHECK(toRead.size() == 2);
            BOOST_CHECK(toWaitFor.size() == 0);
            v1.dismiss();
            mgr.endVisit(visitId - 1, false);
        }
        v1.dismiss();
        v2.dismiss();
    }
    mgr.endVisit(numVisits - 1, false);
}


BOOST_AUTO_TEST_CASE(overlappingVisitsTest) {

    BOOST_TEST_MESSAGE("    - ChunkManager test: sequence of overlapping visits");
    SharedSimpleObjectChunkManager mgr("test");
    // unlink the shared memory object immediately (it remains available until the test process exits)
    SharedSimpleObjectChunkManager::destroyInstance("test");

    // Process a series of overlapping visits
    static int64_t const numVisits = 50;
    for (int64_t visitId = 1; visitId < numVisits; ++visitId) {

        ScopeGuard v1(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, visitId - 1, false));
        ScopeGuard v2(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, visitId, false));
        std::vector<SObjChunk> toRead;
        std::vector<SObjChunk> toWaitFor;
        std::vector<int64_t>   chunkIds;

        mgr.registerVisit(visitId);
        chunkIds.push_back(visitId);
        chunkIds.push_back(visitId + 1);
        mgr.startVisit(toRead, toWaitFor, visitId, chunkIds);
        if (visitId > 1) {
            BOOST_CHECK(toRead.size() == 1 && toRead[0].getId() == visitId + 1);
            BOOST_CHECK(toWaitFor.size() == 1 && toWaitFor[0].getId() == visitId);
            // make sure that waiting for ownership times out (there is chunk overlap,
            // and visitId - 1 hasn't yet been ended)
            TimeSpec deadline;
            deadline.systemTime();
            deadline += 0.02;
            BOOST_CHECK_THROW(mgr.waitForOwnership(toRead, toWaitFor, visitId, deadline), Timeout);
            // end the previous visit
            v1.dismiss();
            mgr.endVisit(visitId - 1, false);
            // now check to make sure visitId can proceed
            mgr.waitForOwnership(toRead, toWaitFor, visitId, deadline);
            // note - since we never actually read() chunks, all chunks are permanently in
            // the unusable state, i.e. they will always be inserted into toRead once ownership
            // is acquired.
            BOOST_CHECK(toRead.size() == 1 && toRead[0].getId() == visitId);
            BOOST_CHECK(toWaitFor.size() == 0);
        } else {
            BOOST_CHECK(toRead.size() == 2);
            BOOST_CHECK(toWaitFor.size() == 0);
            v1.dismiss();
        }
        v2.dismiss();
    }
    mgr.endVisit(numVisits - 1, false);
}

