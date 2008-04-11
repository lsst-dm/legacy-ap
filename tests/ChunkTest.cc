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
#include <string>
#include <vector>

#include <boost/version.hpp>
#include <boost/bind.hpp>
#include <boost/shared_array.hpp>
#if BOOST_VERSION < 103400
#   include <boost/test/auto_unit_test.hpp>
#   define BOOST_TEST_MESSAGE BOOST_MESSAGE
#else
#   include <boost/test/unit_test.hpp>
#endif

#include <lsst/ap/Common.h>
#include <lsst/ap/Chunk.h>
#include <lsst/ap/ChunkManager.h>
#include <lsst/ap/Random.h>
#include <lsst/ap/ScopeGuard.h>

using namespace lsst::ap;


namespace {

typedef SharedSimpleObjectChunkManager::SimpleObjectChunk SObjChunk;


static std::string const makeTempFile() {
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
static SObjChunk const createChunk() {
    SharedSimpleObjectChunkManager mgr("test");
    // unlink the shared memory object immediately (it remains available until the test process exits)
    SharedSimpleObjectChunkManager::destroyInstance("test");

    std::vector<int64_t>   chunkIds;
    std::vector<SObjChunk> toWaitFor;
    std::vector<SObjChunk> toRead;

    initRandom();
    mgr.registerVisit(1);
    int64_t chunkId = static_cast<int64_t>(uniformRandom()*1.0e9);
    chunkIds.push_back(chunkId);
    mgr.startVisit(toRead, toWaitFor, 1, chunkIds);
    BOOST_CHECK_MESSAGE(toWaitFor.size() == 0, "Chunk was already registered with the manager");
    BOOST_CHECK_MESSAGE(toRead.size() == 1 && toRead[0].getId() == chunkId, "Couldn't locate chunk via manager");
    return toRead[0];
}


// Pick num ids at random, with values between min (inclusive) and max (exclusive).
// Then sort the picks and remove duplicates.
static void pickIds(std::vector<int64_t> & ids, size_t const num, int64_t const min, int64_t const max) {
    BOOST_REQUIRE(max > min && min >= 0);
    for (size_t i = 0; i < num; ++i) {
        ids.push_back(min + static_cast<int64_t>(uniformRandom()*(max - min)));
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}


static void appendObjects(SObjChunk & chunk, size_t const num) {
    SimpleObject obj;
    int64_t      id = static_cast<int64_t>(chunk.size());
    for (size_t i = 0; i < num; ++i, ++id) {
        obj._objectId = id;
        // don't care where objects are for this test - fill everything except id with random data
        Point p   = Point::random();
        obj._ra   = p._ra;
        obj._decl = p._dec;
        for (int i = 0; i < lsst::afw::image::Filter::NUM_FILTERS; ++i) {
            obj._varProb[i] = static_cast<int8_t>(uniformRandom()*100.0);
        }
        chunk.insert(obj);
    }
}


static boost::shared_array<SimpleObject> const copyData(SObjChunk & chunk) {
    uint32_t sz = chunk.size();
    boost::shared_array<SimpleObject> copy(new SimpleObject[sz]);
    std::memset(copy.get(), 0, sz*sizeof(SimpleObject));
    for (uint32_t i = 0; i < sz; ++ i) {
        copy[i] = chunk.get(i);
    }
    return copy;
}


static bool isInDelta(SObjChunk const & chunk, uint32_t const i, uint32_t const off = 0) {
    return (chunk.getFlag(i - off) & SObjChunk::IN_DELTA) != 0;
}

static bool isUncommitted(SObjChunk const & chunk, uint32_t const i, uint32_t const off = 0) {
    return (chunk.getFlag(i - off) & SObjChunk::UNCOMMITTED) != 0;
}

static bool isInserted(SObjChunk const & chunk, uint32_t const i, uint32_t const off = 0) {
    return (chunk.getFlag(i - off) & SObjChunk::INSERTED) != 0;
}

static bool isDeleted(SObjChunk const & chunk, uint32_t const i, uint32_t const off = 0) {
    return (chunk.getFlag(i - off) & SObjChunk::DELETED) != 0;
}


static void verifyChunk(
    SObjChunk            const & chunk,
    std::vector<int64_t> const & deletes,
    bool                 const   packed,
    bool                 const   committed
) {
    uint32_t size  = chunk.size();
    uint32_t delta = chunk.delta();
    if (packed) {
        if (delta < size) {
            delta += std::lower_bound(deletes.begin(), deletes.end(), chunk.get(delta).getId()) - deletes.begin();
        } else {
            delta += deletes.size();
        }
        size += deletes.size();
    }

    std::vector<int64_t>::const_iterator gap(deletes.begin());
    int64_t  nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int64_t>::max();
    uint32_t off     = 0;
    bool     foundInsert = false;

    // loop over normal entries
    for (uint32_t i = 0; i < delta; ++i) {
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
            nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int64_t>::max();
        } else {
            BOOST_CHECK_MESSAGE(!isInDelta(chunk, i, off), "entry marked IN_DELTA is below the chunk delta index");
            BOOST_CHECK_MESSAGE(!isInserted(chunk, i, off), "INSERTED entry not marked IN_DELTA");
            BOOST_CHECK_MESSAGE(!isDeleted(chunk, i, off), "entry incorrectly marked as DELETED");
            BOOST_CHECK_MESSAGE(!isUncommitted(chunk, i, off), "uncommitted entry not IN_DELTA must be DELETED");
        }
    }

    // then loop over entries in the chunk delta
    for (uint32_t i = delta; i < size; ++i) {
        if (i == nextGap) {
            if (packed) {
                if (i - off < chunk.size()) {
                    BOOST_CHECK_MESSAGE(chunk.get(i - off).getId() > i, "deleted entry wasn't packed");
                }
                ++off;
                nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int64_t>::max();
                continue;
            }
            BOOST_CHECK_MESSAGE(isDeleted(chunk, i, off), "deleted entry wasn't marked DELETED");
            nextGap = (gap != deletes.end()) ? *gap++ : std::numeric_limits<int64_t>::max();
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


static void verifyData(SObjChunk const & chunk, boost::shared_array<SimpleObject> const & data) {
    uint32_t size = chunk.size();
    for (uint32_t i = 0, b = 0; i < size; ++b) {
        uint32_t n = chunk.entries(b);
        int cmp = std::memcmp(chunk.getBlock(b), &(data[i]), sizeof(SimpleObject)*n);
        i += n;
        BOOST_CHECK_MESSAGE(cmp ==  0, "chunk IO resulted in data corruption");
    }
}


static void verifyPackedData(
    SObjChunk                         const & chunk,
    boost::shared_array<SimpleObject> const & data,
    std::vector<int64_t>              const & deletes
) {
    uint32_t size = chunk.size();
    std::vector<int64_t>::const_iterator gap(deletes.begin());
    for (uint32_t i = 0, j = 0; i < size; ++j) {
        if (gap != deletes.end() && *gap == j) {
            ++gap;
            continue;
        }
        BOOST_CHECK_MESSAGE(data[j] == chunk.get(i), "chunk IO resulted in data corruption: " <<
            data[j]._objectId << " != " << chunk.get(i)._objectId);
        ++i;
    }
}


static void wrCycle(SObjChunk & c, std::vector<int64_t> & d, std::string const & name, bool const compress) {
    uint32_t const size  = c.size();
    uint32_t const delta = c.delta();

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


static void wrdCycle(
    SObjChunk            & c,
    std::vector<int64_t> & d,
    std::string    const & name,
    bool           const   packed,
    bool           const   compress
) {
    uint32_t const size  = c.size();
    uint32_t const delta = c.delta();
    std::string    deltaName(makeTempFile());
    ScopeGuard     guard(boost::bind(::unlink, deltaName.c_str()));

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
    SharedSimpleObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int64_t> v;
    SObjChunk c(createChunk());
    BOOST_CHECK_MESSAGE(c.size() == 0,  "freshly created chunk isn't empty");
    appendObjects(c, static_cast<size_t>(uniformRandom()*32768.0) + 1024);
    BOOST_CHECK_MESSAGE(c.delta() == 0, "freshly created chunk delta should be 0");
    verifyChunk(c, v, false, false);

    boost::shared_array<SimpleObject> data(copyData(c));
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));
    wrCycle(c, v, name, false);
    verifyData(c, data);
    wrCycle(c, v, name, true);
    verifyData(c, data);
}


BOOST_AUTO_TEST_CASE(chunkIoDeltaTest) {
    BOOST_TEST_MESSAGE("    - Chunk IO test (with delta and no deletes)");
    SharedSimpleObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int64_t> v;
    SObjChunk c(createChunk());
    appendObjects(c, static_cast<size_t>(uniformRandom()*16384.0) + 1024);
    verifyChunk(c, v, false, false);

    boost::shared_array<SimpleObject> data(copyData(c));
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));

    wrCycle(c, v, name, false);
    verifyData(c, data);
    uint32_t size = c.size();
    appendObjects(c, static_cast<size_t>(uniformRandom()*8192.0));
    c.rollback();
    BOOST_CHECK(c.size()  == size);
    BOOST_CHECK(c.delta() == size);
    appendObjects(c, static_cast<size_t>(uniformRandom()*8192.0));
    data  = copyData(c);
    wrdCycle(c, v, name, false, false);
    verifyData(c, data);

    // repeat test, this time with compression
    wrCycle(c, v, name, true);
    verifyData(c, data);
    appendObjects(c, static_cast<size_t>(uniformRandom()*8192.0));
    data  = copyData(c);
    wrdCycle(c, v, name, false, true);
    verifyData(c, data);
}


BOOST_AUTO_TEST_CASE(chunkIoDeltaDelTest) {
    BOOST_TEST_MESSAGE("    - Chunk IO test (with delta and deletes)");
    SharedSimpleObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int64_t> v;
    SObjChunk c(createChunk());
    appendObjects(c, static_cast<size_t>(uniformRandom()*65536.0) + 1024);
    verifyChunk(c, v, false, false);

    boost::shared_array<SimpleObject> data(copyData(c));
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));

    wrCycle(c, v, name, false);
    uint32_t size = c.size();
    verifyData(c, data);
    appendObjects(c, static_cast<size_t>(uniformRandom()*8192.0));
    data  = copyData(c);
    uint32_t deltaSize = c.size() - size;
    wrdCycle(c, v, name, false, false);
    verifyData(c, data);

    // pick ids (== index in chunk) of objects to delete
    pickIds(v, static_cast<size_t>(static_cast<double>(size)*0.25*uniformRandom()), 0, size);
    pickIds(v, static_cast<size_t>(static_cast<double>(deltaSize)*0.25*uniformRandom()), size, size + deltaSize);

    // and delete those objects
    for (std::vector<int64_t>::const_iterator i = v.begin(); i != v.end(); ++i) {
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
    SharedSimpleObjectChunkManager mgr("test");
    ScopeGuard forMgr(boost::bind(&SharedSimpleObjectChunkManager::endVisit, &mgr, 1, true));

    std::vector<int64_t> v;
    SObjChunk c(createChunk());
    verifyChunk(c, v, false, false);

    boost::shared_array<SimpleObject> data;
    std::string name(makeTempFile());
    ScopeGuard  guard(boost::bind(::unlink, name.c_str()));

    // test reading/writing of empty chunk and chunk delta files
    wrCycle(c, v, name, false);
    verifyData(c, data);
    wrdCycle(c, v, name, false, false);
    verifyData(c, data);
}

