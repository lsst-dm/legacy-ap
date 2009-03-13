// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Testing of Chunk class functionality (I/O, etc...).
 *
 * @ingroup associate
 */

#include <unistd.h>

#include <cstring>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "boost/bind.hpp"
#include "boost/shared_array.hpp"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ChunkTest
#include "boost/test/unit_test.hpp"

#include "lsst/afw/math/Random.h"

#include "lsst/ap/Common.h"
#include "lsst/ap/Chunk.h"
#include "lsst/ap/ChunkManager.h"
#include "lsst/ap/Point.h"
#include "lsst/ap/ScopeGuard.h"
#include "lsst/ap/Time.h"


using lsst::afw::math::Random;
using namespace lsst::ap;


typedef SharedObjectChunkManager::ObjectChunk ObjChunk;

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


bool coinToss(double p) {
    return rng().uniform() <= p;
}


std::string const makeTempFile() {
    char name[32];
    std::strncpy(name, "/tmp/ChunkTest.XXXXXX", 31);
    name[31] = 0;
    int const fd = ::mkstemp(name);
    if (fd < 1) {
        BOOST_FAIL("Failed to create temporary file for testing purposes");
    }
    ::close(fd);
    return std::string(name);
}


// Creates a single, initially empty, chunk belonging to visit 1.
ObjChunk const createChunk() {
    SharedObjectChunkManager mgr("test");
    // unlink the shared memory object immediately (it remains available until the test process exits)
    SharedObjectChunkManager::destroyInstance("test");

    std::vector<int>      chunkIds;
    std::vector<ObjChunk> toWaitFor;
    std::vector<ObjChunk> toRead;

    mgr.registerVisit(1);
    int chunkId = static_cast<int>(rng().uniformInt(1000000000));
    chunkIds.push_back(chunkId);
    mgr.startVisit(toRead, toWaitFor, 1, chunkIds);
    BOOST_CHECK_MESSAGE(toWaitFor.size() == 0, "Chunk was already registered with the manager");
    BOOST_CHECK_MESSAGE(toRead.size() == 1 && toRead[0].getId() == chunkId, "Couldn't locate chunk via manager");
    return toRead[0];
}


// Pick num ids at random, with values between min (inclusive) and max (exclusive).
// Then sort the picks and remove duplicates.
void pickIds(std::vector<int> & ids, int const num, int const min, int const max) {
    BOOST_REQUIRE(max > min && min >= 0);
    for (int i = 0; i < num; ++i) {
        ids.push_back(static_cast<int>(rng().flat(min, max)));
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}


void appendObjects(ObjChunk & chunk, int const num) {
    Object obj;
    int id = chunk.size();
    for (int i = 0; i < num; ++i, ++id) {
        obj._objectId = id;
        // don't care where objects are for this test - fill everything except id with random data
        Point p     = Point::random(rng());
        obj._ra     = p._ra;
        obj._decl   = p._dec;
        obj._pmRa   = rng().gaussian()*10.0;
        obj._pmDecl = rng().gaussian()*10.0;
        for (int i = 0; i < lsst::afw::image::Filter::NUM_FILTERS; ++i) {
            obj._varProb[i] = static_cast<boost::int16_t>(rng().uniformInt(100));
        }
        chunk.insert(obj);
    }
}


boost::shared_array<Object> const copyData(ObjChunk & chunk) {
    int sz = chunk.size();
    boost::shared_array<Object> copy(new Object[sz]);
    std::memset(copy.get(), 0, sz*sizeof(Object));
    for (int i = 0; i < sz; ++ i) {
        copy[i] = chunk.get(i);
    }
    return copy;
}


bool isInDelta(ObjChunk const & chunk, int const i, int const off = 0) {
    return (chunk.getFlag(i - off) & ObjChunk::IN_DELTA) != 0;
}

bool isUncommitted(ObjChunk const & chunk, int const i, int const off = 0) {
    return (chunk.getFlag(i - off) & ObjChunk::UNCOMMITTED) != 0;
}

bool isInserted(ObjChunk const & chunk, int const i, int const off = 0) {
    return (chunk.getFlag(i - off) & ObjChunk::INSERTED) != 0;
}

bool isDeleted(ObjChunk const & chunk, int const i, int const off = 0) {
    return (chunk.getFlag(i - off) & ObjChunk::DELETED) != 0;
}


void verifyChunk(
    ObjChunk         const & chunk,
    std::vector<int> const & deletes,
    bool             const   packed,
    bool             const   committed
) {
    int size  = chunk.size();
    int delta = chunk.delta();
    if (packed) {
        if (delta < size) {
            delta += std::lower_bound(deletes.begin(), deletes.end(), chunk.get(delta).getId()) - deletes.begin();
        } else {
            delta += deletes.size();
        }
        size += deletes.size();
    }

    std::vector<int>::const_iterator gap(deletes.begin());
    int nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int>::max();
    int off = 0;
    bool foundInsert = false;

    // loop over normal entries
    for (int i = 0; i < delta; ++i) {
        if (i == nextGap) {
            if (packed) {
                if (i - off < chunk.size()) {
                    BOOST_CHECK_MESSAGE(chunk.get(i - off).getId() > i, "deleted entry wasn't packed");
                }
                ++off;
            } else {
                BOOST_CHECK_MESSAGE(!isInDelta(chunk, i, off), "entry marked IN_DELTA is below the chunk delta index");
                BOOST_CHECK_MESSAGE(!isInserted(chunk, i, off), "found INSERTED entry not marked IN_DELTA");
                BOOST_CHECK_MESSAGE(isDeleted(chunk, i, off), "deleted entry " << i << " wasn't marked DELETED");
                if (committed) {
                    BOOST_CHECK_MESSAGE(!isUncommitted(chunk, i, off), "found uncommitted entry");
                }
            }
            nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int>::max();
        } else {
            BOOST_CHECK_MESSAGE(!isInDelta(chunk, i, off), "entry marked IN_DELTA is below the chunk delta index");
            BOOST_CHECK_MESSAGE(!isInserted(chunk, i, off), "INSERTED entry not marked IN_DELTA");
            BOOST_CHECK_MESSAGE(!isDeleted(chunk, i, off), "entry incorrectly marked as DELETED");
            BOOST_CHECK_MESSAGE(!isUncommitted(chunk, i, off), "uncommitted entry not IN_DELTA must be DELETED");
        }
    }

    // then loop over entries in the chunk delta
    for (int i = delta; i < size; ++i) {
        if (i == nextGap) {
            if (packed) {
                if (i - off < chunk.size()) {
                    BOOST_CHECK_MESSAGE(chunk.get(i - off).getId() > i, "deleted entry wasn't packed");
                }
                ++off;
                nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int>::max();
                continue;
            }
            BOOST_CHECK_MESSAGE(isDeleted(chunk, i, off), "deleted entry wasn't marked DELETED");
            nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int>::max();
        } else {
            BOOST_CHECK_MESSAGE(!isDeleted(chunk, i, off), "entry incorrectly marked as DELETED");
        }
        BOOST_CHECK_MESSAGE(isInDelta(chunk, i, off), "found entry in chunk delta not marked IN_DELTA");
        if (committed) {
            BOOST_CHECK_MESSAGE(!isUncommitted(chunk, i, off), "found uncommitted entry");
            BOOST_CHECK_MESSAGE(!isInserted(chunk, i, off), "found entry marked INSERTED in committed chunk");
        } else if (foundInsert) {
            BOOST_CHECK_MESSAGE(isInserted(chunk, i, off), "inserts must be consecutive (at the end of a chunk)");
            BOOST_CHECK_MESSAGE(isUncommitted(chunk, i, off), "found entry marked as INSERTED but not UNCOMMITTED");
        } else if ((foundInsert = isInserted(chunk, i, off))) {
            BOOST_CHECK_MESSAGE(isUncommitted(chunk, i, off), "found entry marked as INSERTED but not UNCOMMITTED");
        }
    }
}


void verifyData(ObjChunk const & chunk, boost::shared_array<Object> const & data) {
    int size = chunk.size();
    for (int i = 0, b = 0; i < size; ++b) {
        int n = chunk.entries(b);
        int cmp = std::memcmp(chunk.getBlock(b), &(data[i]), sizeof(Object)*n);
        i += n;
        BOOST_CHECK_MESSAGE(cmp ==  0, "chunk IO resulted in data corruption");
    }
}


void verifyPackedData(
    ObjChunk                    const & chunk,
    boost::shared_array<Object> const & data,
    std::vector<int>            const & deletes
) {
    int size = chunk.size();
    std::vector<int>::const_iterator gap(deletes.begin());
    for (int i = 0, j = 0; i < size; ++j) {
        if (gap != deletes.end() && *gap == j) {
            ++gap;
            continue;
        }
        BOOST_CHECK_MESSAGE(data[j] == chunk.get(i), "chunk IO resulted in data corruption: " <<
            data[j]._objectId << " != " << chunk.get(i)._objectId);
        ++i;
    }
}


void wrCycle(ObjChunk & c, std::vector<int> & d, std::string const & name, bool const compress) {
    int const size  = c.size();
    int const delta = c.delta();

    c.write(name, true, compress);
    if (coinToss(0.5)) {
        c.commit(false);
        BOOST_CHECK_MESSAGE(c.delta() == delta, "unexpected chunk delta index modification");
    } else {
        c.commit(true);
        BOOST_CHECK_MESSAGE(c.delta() == size, "chunk delta index must equal chunk size following commit(true)");
    }
    verifyChunk(c, d, false, true);
    c.clear();
    c.read(name, compress);
    BOOST_CHECK_MESSAGE(c.delta() == size, "chunk delta index must equal chunk size following chunk read");
    BOOST_CHECK_MESSAGE(c.size() == size, "chunk IO resulted in data corruption");
    verifyChunk(c, d, false, true);
}


void wrdCycle(
    ObjChunk          & c,
    std::vector<int>  & d,
    std::string const & name,
    bool        const   packed,
    bool        const   compress
) {
    int const size  = c.size();
    int const delta = c.delta();
    std::string deltaName(makeTempFile());
    ScopeGuard guard(boost::bind(::unlink, deltaName.c_str()));

    verifyChunk(c, d, packed, false);
    c.writeDelta(deltaName, true, compress);
    if (coinToss(0.5)) {
        c.commit(false);
        BOOST_CHECK_MESSAGE(c.delta() == delta, "unexpected chunk delta index modification");
    } else {
        c.commit(true);
        BOOST_CHECK_MESSAGE(c.delta() == size, "chunk delta index must equal chunk size following commit(true)");
    }
    verifyChunk(c, d, packed, true);
    c.clear();
    c.read(name, compress);
    BOOST_CHECK_MESSAGE(c.size() == delta, "failed to retrieve chunk entries not in delta");
    BOOST_CHECK_MESSAGE(c.delta() == delta, "chunk delta index must equal chunk size following chunk read");
    if (d.size() == 0) {
        verifyChunk(c, d, packed, true);
    }
    c.readDelta(deltaName, compress);
    BOOST_CHECK_MESSAGE(c.size() == size, "chunk IO resulted in data corruption");
    BOOST_CHECK_MESSAGE(c.delta() == delta, "unexpected chunk delta index modification");
    verifyChunk(c, d, packed, true);
}

} // end of anonymous namespace


BOOST_AUTO_TEST_CASE(chunkIoTest) {
    BOOST_TEST_MESSAGE("    - Chunk IO test (without delta)");
    SharedObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int> v;
    ObjChunk c(createChunk());
    BOOST_CHECK_MESSAGE(c.size() == 0,  "freshly created chunk isn't empty");
    appendObjects(c, static_cast<int>(rng().flat(1024, 32768)));
    BOOST_CHECK_MESSAGE(c.delta() == 0, "freshly created chunk delta should be 0");
    verifyChunk(c, v, false, false);

    boost::shared_array<Object> data(copyData(c));
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));
    wrCycle(c, v, name, false);
    verifyData(c, data);
    wrCycle(c, v, name, true);
    verifyData(c, data);
}


BOOST_AUTO_TEST_CASE(chunkIoDeltaTest) {
    BOOST_TEST_MESSAGE("    - Chunk IO test (with delta and no deletes)");
    SharedObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int> v;
    ObjChunk c(createChunk());
    appendObjects(c, static_cast<int>(rng().flat(1024, 16384)));
    verifyChunk(c, v, false, false);

    boost::shared_array<Object> data(copyData(c));
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));

    wrCycle(c, v, name, false);
    verifyData(c, data);
    int size = c.size();
    appendObjects(c, static_cast<int>(rng().flat(0, 8192)));
    c.rollback();
    BOOST_CHECK(c.size()  == size);
    BOOST_CHECK(c.delta() == size);
    appendObjects(c, static_cast<int>(rng().flat(0, 8192)));
    data  = copyData(c);
    wrdCycle(c, v, name, false, false);
    verifyData(c, data);

    // repeat test, this time with compression
    wrCycle(c, v, name, true);
    verifyData(c, data);
    appendObjects(c, static_cast<int>(rng().flat(0, 8192)));
    data  = copyData(c);
    wrdCycle(c, v, name, false, true);
    verifyData(c, data);
}


BOOST_AUTO_TEST_CASE(chunkIoDeltaDelTest) {
    BOOST_TEST_MESSAGE("    - Chunk IO test (with delta and deletes)");
    SharedObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int> v;
    ObjChunk c(createChunk());
    appendObjects(c, static_cast<int>(rng().flat(1024, 65536)));
    verifyChunk(c, v, false, false);

    boost::shared_array<Object> data(copyData(c));
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));

    wrCycle(c, v, name, false);
    int size = c.size();
    verifyData(c, data);
    appendObjects(c, static_cast<int>(rng().flat(0, 8192)));
    data  = copyData(c);
    int deltaSize = c.size() - size;
    wrdCycle(c, v, name, false, false);
    verifyData(c, data);

    // pick ids (== index in chunk) of objects to delete
    pickIds(v, static_cast<int>(static_cast<double>(size)*0.25*rng().uniform()), 0, size);
    pickIds(v, static_cast<int>(static_cast<double>(deltaSize)*0.25*rng().uniform()), size, size + deltaSize);

    // and delete those objects
    for (std::vector<int>::const_iterator i = v.begin(); i != v.end(); ++i) {
        c.remove(*i);
    }
    verifyChunk(c, v, false, false);

    // write things out, with deletes
    wrdCycle(c, v, name, false, false);
    verifyData(c, data);

    // pack the chunk, try again
    c.pack();
    // write out packed version of non-delta entries
    c.write(name, true, false, false);
    wrdCycle(c, v, name, true, false);
    verifyPackedData(c, data, v);
}


BOOST_AUTO_TEST_CASE(emptyChunkTest) {
    BOOST_TEST_MESSAGE("    - Empty chunk IO test");
    SharedObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int> v;
    ObjChunk c(createChunk());
    verifyChunk(c, v, false, false);

    boost::shared_array<Object> data;
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));

    // test reading/writing of empty chunk and chunk delta files
    wrCycle(c, v, name, false);
    verifyData(c, data);
    wrdCycle(c, v, name, false, false);
    verifyData(c, data);
}

