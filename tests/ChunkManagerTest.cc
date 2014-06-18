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
 * @brief   Tests whether ChunkManager handles successive visits
 *          (including overlapping ones) correctly.
 *
 * @ingroup associate
 */

#include <vector>

#include "boost/bind.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ChunkManagerTest
#include "boost/test/unit_test.hpp"

#include "lsst/pex/exceptions.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/Chunk.h"
#include "lsst/ap/ChunkManager.h"
#include "lsst/ap/ScopeGuard.h"
#include "lsst/ap/Time.h"


using namespace lsst::ap;

typedef SharedObjectChunkManager::ObjectChunk ObjChunk;


BOOST_AUTO_TEST_CASE(disjointVisitsTest) {

    BOOST_TEST_MESSAGE("    - ChunkManager test: sequence of disjoint visits");
    SharedObjectChunkManager mgr("test");
    // unlink the shared memory object immediately (it remains available until the test process exits)
    SharedObjectChunkManager::destroyInstance("test");

    // Process a series of non-overlapping visits
    static int const numVisits = 50;
    for (int visitId = 1; visitId < numVisits; ++visitId) {

        ScopeGuard v1(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, visitId - 1, false));
        ScopeGuard v2(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, visitId, false));
        std::vector<ObjChunk> toRead;
        std::vector<ObjChunk> toWaitFor;
        std::vector<int>      chunkIds;

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
    SharedObjectChunkManager mgr("test");
    // unlink the shared memory object immediately (it remains available until the test process exits)
    SharedObjectChunkManager::destroyInstance("test");

    // Process a series of overlapping visits
    static int const numVisits = 50;
    for (int visitId = 1; visitId < numVisits; ++visitId) {

        ScopeGuard v1(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, visitId - 1, false));
        ScopeGuard v2(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, visitId, false));
        std::vector<ObjChunk> toRead;
        std::vector<ObjChunk> toWaitFor;
        std::vector<int>      chunkIds;

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
            BOOST_CHECK_THROW(mgr.waitForOwnership(toRead, toWaitFor, visitId, deadline),
                              lsst::pex::exceptions::TimeoutError);
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

