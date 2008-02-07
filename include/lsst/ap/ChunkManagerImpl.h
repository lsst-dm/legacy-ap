// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Chunk manager helper classes
 *
 * @ingroup associate
 */

#ifndef LSST_AP_CHUNK_MANAGER_IMPL_H
#define LSST_AP_CHUNK_MANAGER_IMPL_H

#include <climits>
#include <iosfwd>

#include <boost/noncopyable.hpp>
#include <boost/static_assert.hpp>
//#include <boost/type_traits/has_nothrow_constructor.hpp>
//#include <boost/type_traits/has_trivial_destructor.hpp>

#include "Common.h"
#include "Bitset.h"
#include "Mutex.h"
#include "Condition.h"
#include "Chunk.h"
#include "DataTraits.h"
#include "Time.h"


namespace lsst {
namespace ap {
namespace detail {

/**
 * @brief  A set of up to @a NumEntries elements of type @a EntryType, hashed by an integer identifier.
 *
 * The hash table implementation is chained and intrusive -- requirements follow:
 *
 * <ul>
 * <li>@a NumEntries must be a positive power of 2</li>
 * <li>@a EntryType must have @code int64_t getId() @endcode and @code void setId(int64_t) @endcode
 *     methods which get/set an int64_t member field and never throw. An id value of -1 is reserved
 *     for invalid entries.</li>
 * <li>@a EntryType must have @code int getNextInChain() @endcode and
 *     @code void setNextInChain(int) @endcode methods which get/set an
 *     int member field and never throw</li>
 * <li>@a EntryType must have a trivial destructor, and a default constructor that doesn't throw</li>
 * </ul>
 */
template <typename EntryType, uint32_t NumEntries>
class HashedSet : private boost::noncopyable {

    BOOST_STATIC_ASSERT(NumEntries > 0 && (NumEntries & (NumEntries - 1)) == 0);
    BOOST_STATIC_ASSERT(NumEntries < INT_MAX);
//    BOOST_STATIC_ASSERT(boost::has_nothrow_constructor<EntryType>::value);
//    BOOST_STATIC_ASSERT(boost::has_trivial_destructor<EntryType>::value);

public :

    HashedSet();

    /// Returns a pointer to the entry with the given identifier, or null if there is no such entry.
    EntryType * find(int64_t const id) { return const_cast<EntryType *>(doFind(id)); }

    /// Returns a pointer to the entry with the given identifier, or null if there is no such entry.
    EntryType const * find(int64_t const id) const { return doFind(id); }

    EntryType * insert(int64_t const id);

    std::pair<EntryType *, bool> findOrInsert(int64_t const id);

    bool erase(int64_t const id);

    /// Returns the number of entries in the set
    uint32_t size() const { return _size; }

    /// Returns the number of additional entries there is space for in the set
    uint32_t space() const { return NumEntries - _size; }

    /**
     * Returns a pointer to the beginning of the underlying array of entries. Not all array entries will
     * correspond to HashedSet entries : invalid entries - those that do not correspond to an entry that
     * has been added to the set via insert() or findOrInsert() - are marked with an id value of -1.
     */
    EntryType       * begin()       { return _entries; }
    EntryType const * begin() const { return _entries; }

    /// Returns a pointer to the end of the underlying array of entries.
    EntryType       * end()       { return _entries + NumEntries; }
    EntryType const * end() const { return _entries + NumEntries; }

private :

    int       _hashTable[2*NumEntries];
    int       _free;
    uint32_t  _size;
    EntryType _entries[NumEntries];

    EntryType const * doFind(int64_t const id) const;
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
template <typename MutexType, typename DataType, typename TraitsType = DataTraits<DataType> >
class BlockAllocator : private boost::noncopyable {

public :

    BlockAllocator(uint8_t const * const ref, size_t const offset);

    size_t allocate();
    void   allocate(size_t       * const blockOffsets, uint32_t const n);
    void   free    (size_t const * const blockOffsets, uint32_t const n);

private :

    typedef Bitset<uint64_t, TraitsType::NUM_BLOCKS> AllocatorType;

    BOOST_STATIC_ASSERT(TraitsType::ENTRIES_PER_BLOCK_LOG2 >= 9);

    static size_t const BLOCK_SIZE =
        (sizeof(DataType) + sizeof(ChunkEntryFlagType)) << TraitsType::ENTRIES_PER_BLOCK_LOG2;

    MutexType     _mutex;
    AllocatorType _allocator;
    size_t const  _offset;
};


/** @brief  State for a single visit to a field of view. */
class LSST_AP_LOCAL Visit {
public :

    Visit() : _id(-1), _next(-1), _failed(false) {}

    bool    failed() const { return _failed; }
    void    setFailed()    { _failed = true; }

    // HashedSet requirements
    int64_t getId         () const           { return _id;   }
    int     getNextInChain() const           { return _next; }
    void    setId         (int64_t const id) { _id   = id;   }
    void    setNextInChain(int     const id) { _next = id;   }

private :

    int64_t _id;
    int     _next;
    bool    _failed;
};


/** @brief  Tracks a set of visits. */
class LSST_AP_LOCAL VisitTracker : public HashedSet<Visit, MAX_VISITS_IN_FLIGHT> {
public :
    bool isValid(int64_t const visitId) const;
    void print(std::ostream & os) const;
    void print(int64_t const visitId, std::ostream & os) const;
};


/**
 * @brief   Helper class for managing chunks of a single type.
 *
 * To be used exclusively by chunk manager implementations.
 */
template <typename MutexType, typename DataType, typename TraitsType = DataTraits<DataType> >
class SubManager : private boost::noncopyable {

public :

    typedef BlockAllocator<MutexType, DataType, TraitsType>   AllocatorType;
    typedef Chunk<AllocatorType, DataType, TraitsType>        ChunkType;
    typedef ChunkDescriptor<TraitsType::MAX_BLOCKS_PER_CHUNK> ChunkDescriptorType;

    // -- fields ----------------

    static uint32_t const NUM_CHUNKS = TraitsType::MAX_CHUNKS_PER_FOV * MAX_VISITS_IN_FLIGHT;

    HashedSet<ChunkDescriptorType, NUM_CHUNKS> _chunks;
    AllocatorType                              _allocator;

    // -- methods ----------------

    SubManager(uint8_t const * const ref, size_t const offset) : _allocator(ref, offset) {}

    /// Returns the number of chunks under management.
    uint32_t size()  const { return _chunks.size();  }
    /// Returns the number of additional chunks that could be handled by this manager.
    uint32_t space() const { return _chunks.space(); }

    void createOrRegisterInterest(
        std::vector<ChunkType>       & toRead,
        std::vector<ChunkType>       & toWaitFor,
        int64_t                const   visitId,
        std::vector<int64_t>   const & chunkIds
    );

    bool checkForOwnership(
        std::vector<ChunkType> & toRead,
        std::vector<ChunkType> & toWaitFor,
        int64_t const            visitId
    );

    void getChunks(
        std::vector<ChunkType>     & chunks,
        std::vector<int64_t> const & chunkIds
    );

    bool relinquishOwnership(
        int64_t      const   visitId,
        bool         const   rollback,
        VisitTracker const & tracker
    );

    void print(std::ostream & os) const;
    void print(int64_t const chunkId, std::ostream & os) const;
    void printVisit(int64_t const visitId, std::ostream & os) const;
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
    typename MutexType,
    typename DataType,
    typename TraitsType = DataTraits<DataType>
>
class ChunkManagerSingleImpl : private boost::noncopyable {

public :

    typedef SubManager<MutexType, DataType, TraitsType> ManagerType;
    typedef typename ManagerType::ChunkType ChunkType;

private :

    /// Returns the offset of the first data block (relative to its manager).
    static size_t blocks() {
        return (sizeof(ChunkManagerSingleImpl) + 511) & ~static_cast<size_t>(511);
    }

public :

    /**
     * Returns the total number of bytes required for a ChunkManagerSingleImpl
     * instance and it's associated pool of memory blocks.
     */
    static size_t size() {
        return blocks() + ChunkType::BLOCK_SIZE * TraitsType::NUM_BLOCKS;
    }

    ChunkManagerSingleImpl();

    bool isVisitInFlight(int64_t const visitId);
    void failVisit      (int64_t const visitId);
    void registerVisit  (int64_t const visitId);

    void startVisit(
        std::vector<ChunkType>       & toRead,
        std::vector<ChunkType>       & toWaitFor,
        int64_t                const   visitId,
        std::vector<int64_t>   const & chunkIds
    );

    void waitForOwnership(
        std::vector<ChunkType> & toRead,
        std::vector<ChunkType> & toWaitFor,
        int64_t          const   visitId,
        TimeSpec         const & deadline
    );

    void getChunks(
        std::vector<ChunkType>     & chunks,
        std::vector<int64_t> const & chunkIds
    );

    bool endVisit(int64_t const visitId, bool const rollback);

    void printVisits(std::ostream & os) const;
    void printChunks(std::ostream & os) const;
    void printVisit(int64_t const visitId, std::ostream & os) const;
    void printChunk(int64_t const chunkId, std::ostream & os) const;

private :

    mutable MutexType    _mutex;
    Condition<MutexType> _ownerCondition;
    VisitTracker         _visits;
    ManagerType          _data;
};


}}} // end of namespace lsst::ap::detail

#endif // LSST_AP_CHUNK_MANAGER_IMPL_H
