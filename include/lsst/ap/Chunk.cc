// -*- lsst-c++ -*-
//
//##====----------------                                ----------------====##/
//
//! \file   Chunk.cc
//
//##====----------------                                ----------------====##/

#ifndef LSST_AP_CHUNK_CC
#define LSST_AP_CHUNK_CC

#include <cstring>
#include <stdexcept>
#include <vector>

#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>

#include <lsst/ap/DataTraits.h>
#include <lsst/ap/Chunk.h>
#include <lsst/ap/io/FileIo.h>


namespace lsst {
namespace ap {


// -- ChunkDescriptor ----------------

template <uint32_t MaxBlocksPerChunk>
void ChunkDescriptor<MaxBlocksPerChunk>::initialize() {

    _chunkId        = -1;
    _visitId        = -1;
    _nextChunk      = -1;
    _usable         = false;
    _numBlocks      = 0;
    _nextBlock      = 0;
    _index          = 0;
    _size           = 0;
    _delta          = 0;
    _curBlockOffset = 0;

    _interestedParties.clear();
    std::memset(_blocks, 0, sizeof(_blocks));
}


// -- Chunk ----------------

/*! Ensures the chunk has space to hold at least \a n entries. */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::reserve(uint32_t const n) {
    if (n > capacity()) {
        if (n > 0x3fffffff) {
            LSST_AP_THROW(LengthError, "requested chunk capacity is too large");
        }
        uint32_t nb = (n + ((1u << ENTRIES_PER_BLOCK_LOG2) - 1)) >> ENTRIES_PER_BLOCK_LOG2;
        uint32_t b  = _descriptor->_numBlocks;
        _allocator->allocate(&(_descriptor->_blocks[b]), nb - b);
        // zero out the newly allocated blocks
        for (; b < nb; ++b) {
            std::memset(
                map(_descriptor->_blocks[b] + (sizeof(ChunkEntryFlagType) << ENTRIES_PER_BLOCK_LOG2)),
                0,
                sizeof(DataType) << ENTRIES_PER_BLOCK_LOG2
            );
        }
        _descriptor->_numBlocks = nb;
    }
}


/*! Inserts the entry into the chunk, allocating memory if necessary. */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::insert(
    DataType           const & data,
    ChunkEntryFlagType const   flags
) {
   uint32_t const block = _descriptor->_nextBlock;
   size_t   off = _descriptor->_curBlockOffset;
   uint32_t i   = _descriptor->_index;

   if (block == 0 || i >= (1u << ENTRIES_PER_BLOCK_LOG2)) {
       // no current block, or current block is full
       if (block >= _descriptor->_numBlocks) {
           if (block >= MAX_BLOCKS) {
               LSST_AP_THROW(LengthError, "maximum number of blocks per chunk exceeded");
           }
           off = _allocator->allocate();
           // zero out the newly allocated block
           std::memset(
               map(off + (sizeof(ChunkEntryFlagType) << ENTRIES_PER_BLOCK_LOG2)),
               0,
               sizeof(DataType) << ENTRIES_PER_BLOCK_LOG2
           );
           _descriptor->_blocks[block] = off;
           _descriptor->_numBlocks     = block + 1;
       } else {
           off = _descriptor->_blocks[block];
       }
       _descriptor->_nextBlock      = block + 1;
       _descriptor->_curBlockOffset = off;
       i = 0;
   }

   // copy in data
   size_t const addr = off + i*sizeof(DataType) + (sizeof(ChunkEntryFlagType) << ENTRIES_PER_BLOCK_LOG2);
   new (reinterpret_cast<DataType *>(map(addr))) DataType(data);

   *reinterpret_cast<ChunkEntryFlagType *>(map(off + i*sizeof(ChunkEntryFlagType))) = flags;
   _descriptor->_index = i + 1;
   ++_descriptor->_size;
}


/*!
    Walks through the chunk beginning at the \a i-th entry and removes all entries
    marked as DELETED. Never throws.

    \param[in] i    The index of the first entry to test for removal.
    \return         \c true if any entries were removed
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
bool Chunk<AllocatorType, DataType, TraitsType>::pack(uint32_t const i) {

    static uint32_t const fb = (sizeof(ChunkEntryFlagType) << ENTRIES_PER_BLOCK_LOG2);

    assert (i < _descriptor->_size);

    uint32_t delta = 0xffffffff;
    uint32_t sz    = _descriptor->_size;
    uint32_t dBlk  = i >> ENTRIES_PER_BLOCK_LOG2;
    size_t   d     = _descriptor->_blocks[dBlk];
    size_t   dEnd  = d + BLOCK_SIZE;
    uint32_t ib    = i & ((1u << ENTRIES_PER_BLOCK_LOG2) - 1);

    ChunkEntryFlagType * __restrict df = map(d + ib*sizeof(ChunkEntryFlagType));
    d += fb + ib*sizeof(DataType);

    uint32_t sBlk = dBlk;
    size_t   s    = d;
    size_t   sEnd = dEnd;
    ChunkEntryFlagType const * __restrict sf = df;

    while (true) {
        for ( ; s < sEnd; s += sizeof(DataType), sf++) {
            if ((*sf & DELETED) != 0) {
                --sz;
                continue;
            }
            if ((*sf & IN_DELTA) != 0 && delta == 0xffffffff) {
                // set delta
                delta = (dBlk << ENTRIES_PER_BLOCK_LOG2) +
                        (d - _descriptor->_blocks[dBlk] - fb)/sizeof(DataType);
            }
            if (df != sf) {
                *df = *sf;
                std::memcpy(map(d), map(s), sizeof(DataType));
            }
            df++;
            d += sizeof(DataType);
            if (d == dEnd) {
                ++dBlk;
                d    = _descriptor->_blocks[dBlk];
                dEnd = d + BLOCK_SIZE;
                df   = map(d);
                d   += fb;
            }
        }
        ++sBlk;
        if (sBlk >= _descriptor->_nextBlock) {
            break;
        }
        s    = _descriptor->_blocks[sBlk];
        sEnd = s + BLOCK_SIZE;
        sf   = map(s);
        s   += fb;
    }

    if (sz < _descriptor->_size) {
        size_t const cb = _descriptor->_blocks[dBlk];
        _descriptor->_curBlockOffset = cb;
        _descriptor->_nextBlock      = dBlk + 1;
        d -= cb + fb;
        _descriptor->_index = d/sizeof(DataType);
        _descriptor->_size  = sz;
        _descriptor->_delta = std::min(sz, delta);
        return true;
    }
    return false;
}


/*!
    Sets flag values for \a n entries in block \a b, starting at the \a i-th entry.

    \param[in] b        The block containing the entries for which flag values are to be set.
    \param[in] flags    The desired flag value of the block entries.
    \param[in] i        The index of the first entry in block \a b for which the flag value is to be set.
    \param[in] n        The number of block entries for which flag values are to be set.
*/
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::setFlags(
    uint32_t           const b,
    ChunkEntryFlagType const flags,
    uint32_t           const i,
    uint32_t           const n
) {
    assert(i + n <= (1u << ENTRIES_PER_BLOCK_LOG2));

    ChunkEntryFlagType * __restrict df = getFlagBlock(b) + i;
    std::memset(df, flags, n);
}


/*!
    Walks through the chunk and undoes any uncommitted inserts or deletes. Never throws.

    \return     \c true if there were any modifications to rollback
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
bool Chunk<AllocatorType, DataType, TraitsType>::rollback() {
    bool mod  = false;

    for (uint32_t b = 0; b < _descriptor->_nextBlock; ++b) {

        size_t const off = _descriptor->_blocks[b];
        ChunkEntryFlagType * const __restrict flags = map(off);
        uint32_t const e = entries(b);

        for (uint32_t i = 0; i < e; ++i) {
            ChunkEntryFlagType f = flags[i];
            if ((f & INSERTED) != 0) {
                // remove all newly inserted entries
                _descriptor->_nextBlock      = b + 1;
                _descriptor->_curBlockOffset = off;
                _descriptor->_index          = i;
                _descriptor->_size           = i + (b << ENTRIES_PER_BLOCK_LOG2);
                return true;
            } else if ((f & (UNCOMMITTED | DELETED)) == (UNCOMMITTED | DELETED)) {
                // undo uncommitted deletes
                flags[i] = f & ~(UNCOMMITTED | DELETED);
                mod = true;
            }
        }
    }
    return mod;
}


/*!
    Marks all uncommitted deletes/inserts as committed. Never throws.

    \param[in] clearDelta   If set to \c true, then the IN_DELTA flag bit is
                            cleared for each entry.
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::commit(bool clearDelta) {

    uint32_t mask = UNCOMMITTED | INSERTED;
    if (clearDelta) {
        mask |= IN_DELTA;
    }
    mask = ~mask;

    for (uint32_t b = 0; b < _descriptor->_nextBlock; ++b) {
        size_t const off = _descriptor->_blocks[b];
        ChunkEntryFlagType * const __restrict flags = map(off);
        uint32_t const e = entries(b);

        for (uint32_t i = 0; i < e; ++i) {
            flags[i] &= mask;
        }
    }

    if (clearDelta) {
        _descriptor->_delta = _descriptor->_size;
    }
}


/*!
    Marks the entries specified by the indexes in the given array as deleted.

    \param[in] deletes      An array containing the indexes of entries to delete.
                            Must be of length at least \a numDeletes.
    \param[in] numDeletes   The number of entries to mark as deleted.
    \param[in] end          Valid entry indexes must be less than this value.
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::applyDeletes(
    uint32_t const * const deletes,
    uint32_t const         numDeletes,
    uint32_t const         end
) {
    // check that all delete indexes are within the specified range
    for (uint32_t i = 0; i < numDeletes; ++i) {
        if (deletes[i] >= end) {
            LSST_AP_THROW(
                IoError,
                "FATAL: binary chunk delta file contains an invalid delete marker - delta NOT applied"
            );
        }
    }

    // ok - apply deletes.
    for (uint32_t i = 0; i < numDeletes; ++i) {
        uint32_t d = deletes[i];
        ChunkEntryFlagType * f = reinterpret_cast<ChunkEntryFlagType *>(map(
            _descriptor->_blocks[d >> ENTRIES_PER_BLOCK_LOG2] +
            (d & ((1u << ENTRIES_PER_BLOCK_LOG2) - 1))*sizeof(ChunkEntryFlagType)
        ));
        *f |= DELETED;
    }
}


static void doRead(io::SequentialReader & reader, uint8_t * dst, size_t dstlen) {
     while (dstlen > 0) {
        size_t nb = reader.read(dst, dstlen);
        assert(nb <= dstlen);
        if (nb == 0) {
            LSST_AP_THROW(IoError, "Unexpected end of file");
        }
        dst    += nb;
        dstlen -= nb;
    }
}


/*!
    Reads the data from the binary chunk file \a name into this chunk. Note that this chunk is
    emptied immediately on entering the function.

    \param name         The name of binary chunk file to read into memory
    \param compressed   Is the binary chunk file compressed (zlib or gzip compression is supported)?
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::read(
    std::string const & name,
    bool        const   compressed
) {
    clear();
    boost::scoped_ptr<io::SequentialReader> reader;
    if (compressed) {
        reader.reset(new io::CompressedFileReader(name));
    } else {
        reader.reset(new io::SequentialFileReader(name));
    }

    if (reader->finished()) {
        // A non-existant chunk file is equivalent to an empty chunk file
        return;
    }

    // read in the header
    BinChunkHeader header;
    doRead(*reader, reinterpret_cast<uint8_t *>(&header), sizeof(BinChunkHeader));

    // check for fishy smells ...
    if (!header.isValid() || header._numDeletes != 0 || header._recordSize != sizeof(DataType)) {
        // lo! blue-fin tuna
        LSST_AP_THROW(
            IoError,
            "Binary chunk file failed sanity checks - file is of the wrong format, "
            "of the wrong record type, or corrupted."
        );
    }

    uint32_t nr = header._numRecords;
    if (nr == 0) {
        return; // nothing to read in
    }
    reserve(nr);

    int32_t  b  = 0;
    uint32_t nd;

    // read in one memory block at a time
    do {
        nd  = std::min((1u << ENTRIES_PER_BLOCK_LOG2), nr);
        nr -= nd;
        doRead(*reader, reinterpret_cast<uint8_t *>(getBlock(b)), nd*sizeof(DataType));
        setFlags(b, 0, 0, nd);
        ++b;
    } while (nr > 0);

    // update chunk size
    _descriptor->_nextBlock      = b;
    _descriptor->_curBlockOffset = _descriptor->_blocks[b - 1];
    _descriptor->_index          = nd;
    uint32_t sz                  = nd + ((b - 1) << ENTRIES_PER_BLOCK_LOG2);
    _descriptor->_size           = sz;
    _descriptor->_delta          = sz;
}


/*!
    Reads the given binary chunk delta file, performing any indicated deletes and appending new
    records to the end of this chunk. This function provides the strong exception safety
    guarantee.

    \param name         The name of binary chunk file to read into memory
    \param compressed   Is the binary chunk file compressed (zlib or gzip compression is supported)?
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::readDelta(
    std::string const & name,
    bool        const   compressed
) {
    boost::scoped_ptr<io::SequentialReader> reader;
    if (compressed) {
        reader.reset(new io::CompressedFileReader(name));
    } else {
        reader.reset(new io::SequentialFileReader(name));
    }

    if (reader->finished()) {
        // A non-existant chunk delta file is equivalent to an empty chunk delta file
        return;
    }

    // read in the header
    BinChunkHeader header;
    doRead(*reader, reinterpret_cast<uint8_t *>(&header), sizeof(BinChunkHeader));

    if (!header.isValid() || header._recordSize != sizeof(DataType)) {
        LSST_AP_THROW(
            IoError,
            "Binary chunk delta file failed sanity checks - file is of the wrong format, "
            "of the wrong record type, or corrupted."
        );
    }

    // read in indexes of records to delete
    boost::scoped_array<uint32_t> deletes(header._numDeletes > 0 ? new uint32_t[header._numDeletes] : 0);
    doRead(*reader, reinterpret_cast<uint8_t *>(deletes.get()), sizeof(uint32_t)*header._numDeletes);

    // read in records to append
    uint32_t nr = header._numRecords;
    if (nr == 0) {
        applyDeletes(deletes.get(), header._numDeletes, _descriptor->_size);
        return; // nothing more to read in
    }
    reserve(nr + _descriptor->_size);

    // fill up the current block (or the first block if there is no current block)
    uint32_t b  = _descriptor->_nextBlock;
    if (b > 0) {
        --b;
    }
    uint32_t i  = _descriptor->_index;
    uint32_t nd = std::min((1u << ENTRIES_PER_BLOCK_LOG2) - i, nr);
    nr -= nd;
    doRead(*reader, reinterpret_cast<uint8_t *>(&getBlock(b)[i]), nd*sizeof(DataType));
    setFlags(b, IN_DELTA, i, nd);
    nd += i;
    ++b;

    // fill remaining blocks
    while (nr > 0) {
        nd  = std::min((1u << ENTRIES_PER_BLOCK_LOG2), nr);
        nr -= nd;
        doRead(*reader, reinterpret_cast<uint8_t *>(getBlock(b)), nd*sizeof(DataType));
        setFlags(b, IN_DELTA, 0, nd);
        ++b;
    }

    uint32_t const sz = ((b - 1) << ENTRIES_PER_BLOCK_LOG2) + nd;
    applyDeletes(deletes.get(), header._numDeletes, sz);

    // update chunk state to reflect additions
    _descriptor->_nextBlock      = b;
    _descriptor->_curBlockOffset = _descriptor->_blocks[b - 1];
    _descriptor->_index          = nd;
    _descriptor->_size           = sz;
}


/*!
    Writes the data from this chunk to a binary file. Note that deleted records will be
    written out as well - to skip them, call pack() immediately before this function. Note
    that uncomitted deletes/inserts and entries flagged as IN_DELTA do not automatically
    have their status cleared - to do this call commit() with \c true as the argument.

    \param name         The name of binary chunk file to write.
    \param overwrite    Should an existing file with the given name be overwritten?
    \param compressed   Should the binary chunk file contents be compressed (a gzip compatible
                        algorithm will be used)?
    \param withDelta    Should entries marked IN_DELTA be written out?
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::write(
    std::string const & name,
    bool        const   overwrite,
    bool        const   compressed,
    bool        const   withDelta
) {
    boost::scoped_ptr<io::SequentialWriter> writer;
    if (compressed) {
        writer.reset(new io::CompressedFileWriter(name, overwrite));
    } else {
        writer.reset(new io::SequentialFileWriter(name, overwrite));
    }

    // create header
    uint32_t nr = withDelta ? size() : delta();
    BinChunkHeader header;
    header._numRecords = nr;
    header._recordSize = sizeof(DataType);

    writer->write(reinterpret_cast<uint8_t const *>(&header), sizeof(BinChunkHeader));
    for (int32_t b = 0; nr > 0; ++b) {
        uint32_t const nd = std::min((1u << ENTRIES_PER_BLOCK_LOG2), nr);
        writer->write(reinterpret_cast<uint8_t const *>(getBlock(b)), nd*sizeof(DataType));
        nr -= nd;
    }
    writer->finish();
}


/*!
    Writes any deletes and inserts in this chunk to a binary delta file named \a name. Note that
    even on successful function return, uncommitted deletes/inserts are \b not marked as committed
    (allowing for writeDelta calls on multiple chunks to be aggregated into a single transaction).
    Finally, note that deleted delta records are written out by this function. To skip them, call
    pack() on this chunk, passing delta() as the parameter.

    \param name         The name of binary chunk file to write
    \param overwrite    Should an existing file with the given name be overwritten?
    \param compressed   Should the binary chunk delta file contents be compressed (a gzip
                        compatible algorithm will be used)?
 */
template <typename AllocatorType, typename DataType, typename TraitsType>
void Chunk<AllocatorType, DataType, TraitsType>::writeDelta(
    std::string const & name,
    bool        const   overwrite,
    bool        const   compressed
) {
    boost::scoped_ptr<io::SequentialWriter> writer;
    if (compressed) {
        writer.reset(new io::CompressedFileWriter(name, overwrite));
    } else {
        writer.reset(new io::SequentialFileWriter(name, overwrite));
    }

    // collect all deletes (there are likely to be very few, if any)
    std::vector<uint32_t> deletes;

    uint32_t const nb = _descriptor->_nextBlock;
    for (uint32_t b = 0; b < nb; ++b) {
        ChunkEntryFlagType const * const __restrict f = getFlagBlock(b);
        uint32_t const nd = (b < nb - 1) ? (1u << ENTRIES_PER_BLOCK_LOG2) : _descriptor->_index;
        for (uint32_t i = 0; i < nd; ++i) {
            if ((f[i] & DELETED) != 0) {
                deletes.push_back(i + (b << ENTRIES_PER_BLOCK_LOG2));
            }
        }
    }

    // create and write header
    uint32_t fd = _descriptor->_delta;
    BinChunkHeader header;
    header._numRecords = _descriptor->_size - fd;
    header._numDeletes = deletes.size();
    header._recordSize = sizeof(DataType);
    writer->write(reinterpret_cast<uint8_t *>(&header), sizeof(BinChunkHeader));

    // write array of delete indexes
    if (!deletes.empty()) {
        writer->write(reinterpret_cast<uint8_t *>(&deletes.front()), deletes.size()*sizeof(uint32_t));
    }

    // write out IN_DELTA chunk entries
    if (fd < _descriptor->_size) {
        uint32_t b = fd >> ENTRIES_PER_BLOCK_LOG2;
        fd &= (1u << ENTRIES_PER_BLOCK_LOG2) - 1;
        uint32_t nd = (b < nb - 1) ? (1u << ENTRIES_PER_BLOCK_LOG2) : _descriptor->_index;
        writer->write(reinterpret_cast<uint8_t *>(&getBlock(b)[fd]), (nd - fd)*sizeof(DataType));
        for (++b; b < nb; ++b) {
            nd = (b < nb - 1) ? (1u << ENTRIES_PER_BLOCK_LOG2) : _descriptor->_index;
            writer->write(reinterpret_cast<uint8_t *>(getBlock(b)), nd*sizeof(DataType));
        }
    }

    // all done
    writer->finish();
}


}}  // end of namespace lsst::ap

#endif // LSST_AP_CHUNK_CC
