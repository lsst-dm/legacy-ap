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
 * @brief   Chunk manager helper classes
 *
 * @ingroup ap
 */

#ifndef LSST_AP_CHUNK_MANAGER_IMPL_H
#define LSST_AP_CHUNK_MANAGER_IMPL_H

#include <climits>
#include <iosfwd>

#include "boost/noncopyable.hpp"
#include "boost/static_assert.hpp"
//#include "boost/type_traits/has_nothrow_constructor.hpp"
//#include "boost/type_traits/has_trivial_destructor.hpp"

#include "Common.h"
#include "Bitset.h"
#include "Mutex.h"
#include "Condition.h"
#include "Chunk.h"
#include "DataTraits.h"
#include "Time.h"


namespace lsst { namespace ap { namespace detail {

/**
 * @brief  A set of up to @a NumEntries elements of type @a EntryT, hashed by an integer identifier.
 *
 * The hash table implementation is chained and intrusive -- requirements follow:
 *
 * <ul>
 * <li>@a NumEntries must be a positive power of 2</li>
 * <li>@a EntryT must have @code int64_t getId() @endcode and @code void setId(int64_t) @endcode
 *     methods which get/set an int64_t member field and never throw. An id value of -1 is reserved
 *     for invalid entries.</li>
 * <li>@a EntryT must have @code int getNextInChain() @endcode and
 *     @code void setNextInChain(int) @endcode methods which get/set an
 *     int member field and never throw</li>
 * <li>@a EntryT must have a trivial destructor, and a default constructor that doesn't throw</li>
 * </ul>
 */
template <typename EntryT, int NumEntries>
class HashedSet : private boost::noncopyable {
    BOOST_STATIC_ASSERT(NumEntries > 0 && (NumEntries & (NumEntries - 1)) == 0);
    BOOST_STATIC_ASSERT(NumEntries < INT_MAX);
//    BOOST_STATIC_ASSERT(boost::has_nothrow_constructor<EntryT>::value);
//    BOOST_STATIC_ASSERT(boost::has_trivial_destructor<EntryT>::value);

public :
    HashedSet();

    /// Returns a pointer to the entry with the given identifier, or null if there is no such entry.
    EntryT * find(int const id) {
        return const_cast<EntryT *>(doFind(id));
    }

    /// Returns a pointer to the entry with the given identifier, or null if there is no such entry.
    EntryT const * find(int const id) const {
        return doFind(id);
    }

    EntryT * insert(int const id);

    std::pair<EntryT *, bool> findOrInsert(int const id);

    bool erase(int const id);

    /** Returns the number of entries in the set */
    int size() const {
        return _size;
    }

    /** Returns the number of additional entries there is space for in the set */
    int space() const {
        return NumEntries - _size;
    }

    /**
     * Returns a pointer to the beginning of the underlying array of entries. Not all array entries will
     * correspond to HashedSet entries : invalid entries - those that do not correspond to an entry that
     * has been added to the set via insert() or findOrInsert() - are marked with an id value of -1.
     */
    EntryT * begin() {
        return _entries;
    }
    EntryT const * begin() const {
        return _entries;
    }

    /** Returns a pointer to the end of the underlying array of entries. */
    EntryT * end() {
        return _entries + NumEntries;
    }
    EntryT const * end() const {
        return _entries + NumEntries;
    }

private :
    int _hashTable[2*NumEntries];
    int _free;
    int _size;
    EntryT _entries[NumEntries];

    EntryT const * doFind(int const id) const;
};


/**
 * @brief  A thread-safe memory block allocator that uses a Bitset to track which blocks (out of a
 *         fixed size pool of blocks) are in-use/free.
 *
 * The allocator never returns raw pointers - instead, blocks are identified by an offset in bytes
 * relative to the address of the allocator instance.
 *
 * This scheme allows an allocator instance, the memory blocks it manages, and offsets referencing
 * them to be stored in shared memory. Clients then map these offsets to an actual block address
 * simply by adding the offsets to the (process-specific) block allocator address.
 */
template <typename MutexT, typename DataT, typename TraitsT = DataTraits<DataT> >
class BlockAllocator : private boost::noncopyable {
public :
    BlockAllocator(unsigned char const * const ref, std::size_t const offset);

    std::size_t allocate();
    void allocate(std::size_t * const blockOffsets, int const n);
    void free(std::size_t const * const blockOffsets, int const n);

private :
    typedef Bitset<boost::uint64_t, TraitsT::NUM_BLOCKS> Allocator;

    BOOST_STATIC_ASSERT(TraitsT::ENTRIES_PER_BLOCK_LOG2 >= 9);

    static std::size_t const BLOCK_SIZE =
        (sizeof(DataT) + sizeof(ChunkEntryFlag)) << TraitsT::ENTRIES_PER_BLOCK_LOG2;

    MutexT _mutex;
    Allocator _allocator;
    std::size_t const _offset;
};


/** @brief  State for a single visit to a field of view. */
class Visit {
public :
    Visit() : _id(-1), _next(-1), _failed(false) {}

    bool failed() const {
        return _failed;
    }
    void setFailed() {
        _failed = true;
    }

    // HashedSet requirements
    int getId() const {
        return _id;
    }
    int getNextInChain() const {
        return _next;
    }
    void setId(int const id) {
        _id = id;
    }
    void setNextInChain(int const id) {
        _next = id;
    }

private :
    int _id;
    int _next;
    bool _failed;
};


/** @brief  Tracks a set of visits. */
class VisitTracker : public HashedSet<Visit, MAX_VISITS_IN_FLIGHT> {
public :
    bool isValid(int const visitId) const;
    void print(std::ostream & os) const;
    void print(int const visitId, std::ostream & os) const;
};


/**
 * @brief   Helper class for managing chunks of a particular type.
 *
 * To be used exclusively by chunk manager implementations.
 */
template <typename MutexT, typename DataT, typename TraitsT = DataTraits<DataT> >
class SubManager : private boost::noncopyable {
public :
    typedef BlockAllocator<MutexT, DataT, TraitsT> Allocator;
    typedef ChunkRef<Allocator, DataT, TraitsT> Chunk;
    typedef ChunkDescriptor<TraitsT::MAX_BLOCKS_PER_CHUNK> Descriptor;

    // -- fields ----------------

    static int const NUM_CHUNKS = TraitsT::MAX_CHUNKS_PER_FOV * MAX_VISITS_IN_FLIGHT;

    HashedSet<Descriptor, NUM_CHUNKS> _chunks;
    Allocator _allocator;

    // -- methods ----------------

    SubManager(unsigned char const * const ref, std::size_t const offset) : _allocator(ref, offset) {}

    /// Returns the number of chunks under management.
    int size()  const { return _chunks.size();  }
    /// Returns the number of additional chunks that could be handled by this manager.
    int space() const { return _chunks.space(); }

    void createOrRegisterInterest(
        std::vector<Chunk> & toRead,
        std::vector<Chunk> & toWaitFor,
        int const visitId,
        std::vector<int> const & chunkIds
    );

    bool checkForOwnership(
        std::vector<Chunk> & toRead,
        std::vector<Chunk> & toWaitFor,
        int const visitId
    );

    void getChunks(
        std::vector<Chunk> & chunks,
        std::vector<int> const & chunkIds
    );

    bool relinquishOwnership(
        int const visitId,
        bool const rollback,
        VisitTracker const & tracker
    );

    void print(std::ostream & os) const;
    void print(int const chunkId, std::ostream & os) const;
    void printVisit(int const visitId, std::ostream & os) const;
};


/**
 * @brief   A manager for a set of chunks of a single type.
 *
 * Each instance is intended to be a header of (meaning: located in the initial bytes of) a
 * large contiguous memory block M, which it then takes charge of managing. In particular, M
 * will contain all chunk descriptors and data, as well as all data structures required for
 * management tasks (visit tracking, memory allocation, and synchronization). M is initialized
 * by calling the placement new operator with the address of M as a parameter. To avoid any
 * data alignment issues, M should begin at an address that is a multiple of 16 bytes (or
 * some larger power of 2).
 *
 * Finally, note that instances of this class do not contain a single pointer - when necessary, offsets
 * (relative to some known address, e.g. of the manager instance itself) are stored instead. This, in
 * conjunction with an appropriate choice of mutex type, makes the class suitable for placement into
 * shared memory.
 */
template <
    typename MutexT,
    typename DataT,
    typename TraitsT = DataTraits<DataT>
>
class ChunkManagerImpl : private boost::noncopyable {
public :
    typedef SubManager<MutexT, DataT, TraitsT> Manager;
    typedef typename Manager::Chunk Chunk;

private :
    /// Returns the offset of the first data block (relative to its manager).
    static std::size_t blocks() {
        return (sizeof(ChunkManagerImpl) + 511) & ~static_cast<size_t>(511);
    }

public :
    /**
     * Returns the total number of bytes required for a ChunkManagerImpl
     * instance and it's associated pool of memory blocks.
     */
    static std::size_t size() {
        return blocks() + Chunk::BLOCK_SIZE * TraitsT::NUM_BLOCKS;
    }

    ChunkManagerImpl();

    bool isVisitInFlight(int const visitId);
    void failVisit(int const visitId);
    void registerVisit(int const visitId);

    void startVisit(
        std::vector<Chunk> & toRead,
        std::vector<Chunk> & toWaitFor,
        int const visitId,
        std::vector<int> const & chunkIds
    );

    void waitForOwnership(
        std::vector<Chunk> & toRead,
        std::vector<Chunk> & toWaitFor,
        int const visitId,
        TimeSpec const & deadline
    );

    void getChunks(
        std::vector<Chunk> & chunks,
        std::vector<int> const & chunkIds
    );

    bool endVisit(int const visitId, bool const rollback);

    void printVisits(std::ostream & os) const;
    void printChunks(std::ostream & os) const;
    void printVisit(int const visitId, std::ostream & os) const;
    void printChunk(int const chunkId, std::ostream & os) const;

private :
    void rollbackAllExcept(int const visitId);

    mutable MutexT    _mutex;
    Condition<MutexT> _ownerCondition;
    VisitTracker      _visits;
    Manager           _data;
};


}}} // end of namespace lsst::ap::detail

#endif // LSST_AP_CHUNK_MANAGER_IMPL_H
