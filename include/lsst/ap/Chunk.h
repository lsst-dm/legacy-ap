// -*- lsst-c++ -*-

/**
 * @file
 * @brief   Classes for holding spatial chunks of things in memory.
 *
 * @ingroup associate
 */

#ifndef LSST_AP_CHUNK_H
#define LSST_AP_CHUNK_H

#include <cassert>
#include <stdexcept>
#include <algorithm>

#include <boost/noncopyable.hpp>
#include <boost/static_assert.hpp>
// #include <boost/type_traits/has_trivial_assign.hpp>
// #include <boost/type_traits/has_trivial_copy.hpp>
// #include <boost/type_traits/has_trivial_destructor.hpp>

#include "Common.h"
#include "DataTraits.h"
#include "Fifo.h"


namespace lsst {
namespace ap {


/** @brief  Simple header for binary chunk files -- allows some sanity checking at read time. */
struct BinChunkHeader {

    static uint32_t const MAGIC = 0xdecade14;

    uint32_t _magic;
    uint32_t _numRecords;
    uint32_t _numDeletes;
    uint32_t _recordSize;

    BinChunkHeader() :
        _magic(MAGIC),
        _numRecords(0),
        _numDeletes(0),
        _recordSize(0)
    {}

    bool isValid() const {
        return _magic == MAGIC;
    }
};


/**
 * @brief  A generic descriptor containing state for different kinds of chunks.
 *
 * State is data and memory type agnostic (that is, the structure contains no
 * pointers and can therefore safely be placed in shared memory).
 */
template <uint32_t MaxBlocksPerChunk>
class ChunkDescriptor : private boost::noncopyable {

public :

    int64_t  _chunkId;   ///< Identifier for the chunk
    int64_t  _visitId;   ///< Identifier for the visit that currently owns the chunk
    int      _nextChunk; ///< Index of the next chunk in the same hash bucket as this one

    /**
     * Flag indicating whether the chunk is in a usable state. A chunk is usable once
     * its data has been successfully and completely read from permanent storage.
     */
    bool     _usable;

    uint32_t _numBlocks;      ///< Number of memory blocks allocated
    uint32_t _nextBlock;      ///< Index of the next block to insert into
    uint32_t _index;          ///< Index of the next free entry in the current block
    uint32_t _size;           ///< Total number of entries
    uint32_t _delta;          ///< Index of first entry marked IN_DELTA
    size_t   _curBlockOffset; ///< Offset of the current block

    /// FIFO of visits to a FOV that overlaps the chunk
    Fifo<MAX_VISITS_IN_FLIGHT> _interestedParties;
    /// List of memory block offsets for allocated blocks
    size_t _blocks[MaxBlocksPerChunk];


    ChunkDescriptor() { initialize(); }

    void initialize();

    void clear() {
        _nextBlock      = 0;
        _index          = 0;
        _size           = 0;
        _delta          = 0;
        _curBlockOffset = 0;
    }

    bool operator<(ChunkDescriptor const & cd) const {
        return _visitId < cd._visitId || (_visitId == cd._visitId && _chunkId < cd._chunkId);
    }

    // detail::HashedSet requirements
    int64_t getId()          const           { return _chunkId;   }
    int     getNextInChain() const           { return _nextChunk; }
    void    setId         (int64_t const id) { _chunkId   = id;   }
    void    setNextInChain(int     const id) { _nextChunk = id;   }
};


/**
 * @def ChunkEntryFlagType
 *
 * The smallest unsigned integral type into which the values of the Chunk::EntryFlag enum
 * can fit. This typedef controls the size of the flag word -- simply using an array of
 * Chunk::EntryFlag values results in sizeof(int) bytes per flag, which is excessive.
 */
typedef uint8_t ChunkEntryFlagType;


/**
 * @brief Provides access to a chunk of things in a rectangular (in right ascension/declination)
 *        region of the sky.
 *
 * This class is essentially a data type specific wrapper for a raw ChunkDescriptor that handles
 * the details of mapping entries and flags to memory blocks provided by a chunk manager. Physical
 * chunk data is referred to via pointers: copy construction and assignment is shallow, and destruction
 * doesn't remove the underlying chunk from memory. Chunk instances can therefore be treated like
 * references.
 *
 * The chunk entry type (a.k.a data type) @b must support a copy constructor that never throws,
 * and @b cannot contain any pointers or references. This is because in some implementations chunk
 * entries are placed in shared memory (which can be mapped to different virtual addresses in
 * different processes).
 *
 * A copy-on-write approach to modifications is enforced: instead of in-place modification,
 * users should mark an existing entry as removed, then insert a copy of a suitably modified entry.
 *
 * <h2>Chunk Memory Layout</h2>
 *
 * A chunk consists of multiple fixed size memory blocks, each of which contains ENTRIES_PER_BLOCK
 * (a power of 2) data items as well as ENTRIES_PER_BLOCK flag words (used to track the state of
 * a chunk entry). A chunk is grown by adding memory blocks without touching existing blocks.
 *
 * This means memory for chunks can be allocated using a trivial fixed size block allocator,
 * and as a consequence:
 * <ul>
 * <li> growing the chunk is O(1) -- allocation never causes a memcpy of the entire chunk </li>
 * <li> freeing a chunk on the other hand is O(N). This isn't a big problem since N is small
 *      (the number of memory blocks belonging to the chunk) </li>
 * <li> the memory allocator can avoid the issue of heap fragmentation entirely. </li>
 * </ul>
 *
 * A single block is laid out as follows:
 * <ul>
 * <li>ENTRIES_PER_BLOCK flag words</li>
 * <li>ENTRIES_PER_BLOCK chunk entries</li>
 * </ul>
 *
 * Note that ENTRIES_PER_BLOCK is a power of 2 greater than or equal to 512. This guarantees that
 * flag words and chunk entries will be located in seperate cache-lines on all modern CPUs (supposing
 * the block itself is properly aligned to a 128byte boundary), and that chunk entries will not suffer
 * from data alignment issues. Also, note that the section of each block dedicated to chunk entries is
 * zeroed on allocation - any padding bytes in a chunk entry should therefore show up as zeroes in
 * memory/when written to disk.
 *
 * Flags for each element of the chunk are stored contiguously and seperately from the
 * actual data objects because:
 *
 * <ul>
 * <li> this allows memory efficient scans of the flag list. When looking for a few data objects (say
 *      those flagged as deleted) among many, this means every memory access retrieves useful data. </li>
 * <li> it seperates internal (non-persistent) bookkeeping data from the data object, which need
 *      only contain fields for things that must be persisted. This means one can sometimes eliminate
 *      some dead space when writing chunk memory to disk. </li>
 * </ul>
 *
 * The following diagram illustrates how the chunk entries are organized and flagged (discontiguous
 * blocks have been visually pasted together to form a logical array of chunk entries):
 *
 * <pre>
 * | block 0| block 1| block 2| block 3| block 4| block 5| block 6| block 7| block 8|
 * |        |        |        |        |        |        |        |        |        |
 * |********|********|********|********|******==|========|========|========|==      |
 * |       -|        |        |-    -  |        |    -   |   -    |        |        |
 * |        |        |        |        |        |        |        |    ++++|++      |
 * |       ~|        |        |        |        |    ~   |        |    ~~~~|~~      |
 * </pre>
 *
 * Legend
 * <dl>
 * <dt> '*' </dt>
 * <dd> Indicates an unmodified chunk entry (since the last time it was read from permanent storage)</dd>
 * <dt> '=' IN_DELTA </dt>
 * <dd> A chunk entry that was inserted during the course of nightly processing. These are the entries
 *      that will be written out to chunk delta files (see below). </dd>
 * <dt> '-' DELETED </dt>
 * <dd> Marks a chunk entry as deleted </dd>
 * <dt> '+' INSERTED </dt>
 * <dd> Marks a chunk entry as inserted </dd>
 * <dt> '~' UNCOMMITTED </dt>
 * <dd> Marks the action associated with the chunk entry (delete or insert) as uncomitted, i.e.
 *      as not yet safely recorded on disk. </dd>
 * </dl>
 *
 * Note the following properties:
 * <ul>
 * <li> All entries marked IN_DELTA ('=') are (logically) contiguous and are to be found
 *      at the end of the (logical) chunk entry array. </li>
 * <li> An entry can be marked INSERTED ('+') only if it is also marked IN_DELTA ('=') </li>
 * <li> All entries marked INSERTED ('+') are (logically) contiguous and are to be found
 *      at the end of the (logical) chunk entry array. </li>
 * <li> All entries marked INSERTED ('+') are also marked UNCOMMITTED ('~') </li>
 * </ul>
 *
 * This strategy allows in-memory rollback of a set of accumulated changes (inserts and deletes),
 * fast location of classes of entries (via memory efficient flag scans), and the ability to record
 * changes to a chunk relative to some initial state. This avoids having to write out entire chunks
 * to record a small number of incremental modifications.
 */
template <typename AllocatorType, typename DataType, typename TraitsType = DataTraits<DataType> >
class Chunk {

//    BOOST_STATIC_ASSERT(boost::has_trivial_assign<DataType>::value);
//    BOOST_STATIC_ASSERT(boost::has_trivial_copy<DataType>::value);
//    BOOST_STATIC_ASSERT(boost::has_trivial_destructor<DataType>::value);

public  :

    typedef DataType EntryType;

    enum EntryFlag {
        IN_DELTA         = 0x01,
        UNCOMMITTED      = 0x02,
        INSERTED         = 0x04,
        DELETED          = 0x08
    };

    static uint32_t const ENTRIES_PER_BLOCK_LOG2 = TraitsType::ENTRIES_PER_BLOCK_LOG2;
    static uint32_t const MAX_BLOCKS             = TraitsType::MAX_BLOCKS_PER_CHUNK;
    static size_t   const BLOCK_SIZE =
        (sizeof(DataType) + sizeof(ChunkEntryFlagType)) << ENTRIES_PER_BLOCK_LOG2;

    typedef ChunkDescriptor<MAX_BLOCKS> DescriptorType;

    Chunk(DescriptorType * desc, AllocatorType * all) : _descriptor(desc), _allocator(all) {}

    ~Chunk() {
        _descriptor = 0;
        _allocator  = 0;
    }

    int64_t getVisitId() const { return _descriptor->_visitId; }
    int64_t getId()      const { return _descriptor->_chunkId; }
    bool    isUsable()   const { return _descriptor->_usable;  }
    void    setUsable()        { _descriptor->_usable = true;  }

    /// Returns a reference to the @a i-th chunk entry.
    DataType const & get(uint32_t const i) const {
        assert(i < _descriptor->_size);
        return *reinterpret_cast<DataType const *>(map(
            _descriptor->_blocks[i >> ENTRIES_PER_BLOCK_LOG2] +
            sizeof(ChunkEntryFlagType)*(1u << ENTRIES_PER_BLOCK_LOG2) +
            (i & ((1u << ENTRIES_PER_BLOCK_LOG2) - 1))*sizeof(DataType)
        ));
    }

    /// Returns a pointer to the @a b-th data block.
    DataType const * getBlock(uint32_t const b) const {
        assert(b < _descriptor->_numBlocks);
        return reinterpret_cast<DataType const *>(map(
            _descriptor->_blocks[b] + (sizeof(ChunkEntryFlagType) << ENTRIES_PER_BLOCK_LOG2)
        ));
    }

    DataType * getBlock(uint32_t const b) {
        assert(b < _descriptor->_numBlocks);
        return reinterpret_cast<DataType *>(map(
            _descriptor->_blocks[b] + (sizeof(ChunkEntryFlagType) << ENTRIES_PER_BLOCK_LOG2)
        ));
    }

    /// Returns the value of the @a i-th flag word.
    ChunkEntryFlagType getFlag(uint32_t const i) const {
        assert(i < _descriptor->_size);
        return *reinterpret_cast<ChunkEntryFlagType const *>(map(
            _descriptor->_blocks[i >> ENTRIES_PER_BLOCK_LOG2] +
            (i & ((1u << ENTRIES_PER_BLOCK_LOG2) - 1))*sizeof(ChunkEntryFlagType)
        ));
    }

    /// Returns a pointer to the @a b-th flag block.
    ChunkEntryFlagType const * getFlagBlock(uint32_t const b) const {
        assert(b < _descriptor->_numBlocks);
        return reinterpret_cast<ChunkEntryFlagType const *>(map(_descriptor->_blocks[b]));
    }

    ChunkEntryFlagType * getFlagBlock(uint32_t const i) {
        assert(i < _descriptor->_numBlocks);
        return reinterpret_cast<ChunkEntryFlagType *>(map(_descriptor->_blocks[i]));
    }

    /// Inserts the given data entry into the chunk, allocating memory if necessary.
    void insert(DataType const & data) { insert(data, IN_DELTA | UNCOMMITTED | INSERTED); }

    /// Marks the i-th entry in the chunk as removed. Never throws.
    void remove(uint32_t const i) {
        assert(i < _descriptor->_size);
        ChunkEntryFlagType * f = reinterpret_cast<ChunkEntryFlagType *>(map(
            _descriptor->_blocks[i >> ENTRIES_PER_BLOCK_LOG2] +
            (i & ((1u << ENTRIES_PER_BLOCK_LOG2) - 1))*sizeof(ChunkEntryFlagType)
        ));
        if ((*f & DELETED) == 0) { // removing an entry that is already deleted is a no-op
            *f |= DELETED | UNCOMMITTED;
        }
    }

    /// Returns the number of entries in the chunk
    uint32_t size() const { return _descriptor->_size; }

    /**
     * Returns the index of the first record in the chunk delta or some value greater than or
     * equal to size() if there are none.
     */
    uint32_t delta() const { return _descriptor->_delta; }

    /// Returns the total number of entries that can be placed in the chunk
    uint32_t capacity() const { return _descriptor->_numBlocks << ENTRIES_PER_BLOCK_LOG2; }

    /// Returns the number of non-empty blocks in the chunk
    uint32_t blocks() const { return _descriptor->_nextBlock; }

    /// Returns the number of entries in the i-th block
    uint32_t entries(uint32_t const i) const {
        assert(i < _descriptor->_nextBlock);
        return (i < _descriptor->_nextBlock - 1) ? (1u << ENTRIES_PER_BLOCK_LOG2) : _descriptor->_index;
    }

    /// Empties the chunk (without deallocating/shrinking memory).
    void clear() { _descriptor->clear(); }

    void reserve(uint32_t const n);

    bool pack(uint32_t const i = 0);

    bool rollback();

    void commit(bool clearDelta = false);

    void read(
        std::string const & name,
        bool        const   compressed
    );

    void readDelta(
        std::string const & name,
        bool        const   compressed
    );

    void write(
        std::string const & name,
        bool        const   overwrite,
        bool        const   compressed,
        bool        const   withDelta = true
    );

    void writeDelta(
        std::string const & name,
        bool        const   overwrite,
        bool        const   compressed
    );

private :

    BOOST_STATIC_ASSERT(TraitsType::ENTRIES_PER_BLOCK_LOG2 >= 9);

    DescriptorType * _descriptor;
    AllocatorType  * _allocator;

    /// Maps the given offset to an actual address.
    uint8_t * map(size_t const off) {
        return reinterpret_cast<uint8_t *>(_allocator) + off;
    }

    /// Maps the given offset to an actual address.
    uint8_t const * map(size_t const off) const {
        return reinterpret_cast<uint8_t const *>(_allocator) + off;
    }

    void applyDeletes(
        uint32_t const * const deletes,
        uint32_t const         numDeletes,
        uint32_t const         end
    );

    void insert(
        DataType           const & data,
        ChunkEntryFlagType const   flags
    );

    void setFlags(
        uint32_t           const b,
        ChunkEntryFlagType const flags,
        uint32_t           const i = 0,
        uint32_t           const n = (1u << ENTRIES_PER_BLOCK_LOG2)
    );
};


}} // end of namespace lsst::ap

#endif // LSST_AP_CHUNK_H

