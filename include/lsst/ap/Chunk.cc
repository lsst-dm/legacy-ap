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
 * @brief   Chunk class implementation.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_CHUNK_CC
#define LSST_AP_CHUNK_CC

#include <cstring>
#include <stdexcept>
#include <vector>

#include "boost/scoped_array.hpp"
#include "boost/scoped_ptr.hpp"

#include "lsst/pex/exceptions.h"

#include "DataTraits.h"
#include "Chunk.h"
#include "io/FileIo.h"


// -- ChunkDescriptor ----------------

template <int MaxBlocksPerChunk>
void lsst::ap::ChunkDescriptor<MaxBlocksPerChunk>::initialize() {

    _chunkId   = -1;
    _visitId   = -1;
    _nextChunk = -1;
    _usable    = false;
    _numBlocks = 0;
    _nextBlock = 0;
    _index     = 0;
    _size      = 0;
    _delta     = 0;
    _curBlockOffset = 0;

    _interestedParties.clear();
    std::memset(_blocks, 0, sizeof(_blocks));
}


// -- ChunkRef ----------------

/** Ensures the chunk has space to hold at least @a n entries. */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::reserve(int const n) {
    if (n > capacity()) {
        if (n > 0x3fffffff) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                "Requested chunk capacity is too large");
        }
        int nb = (n + ((1 << ENTRIES_PER_BLOCK_LOG2) - 1)) >> ENTRIES_PER_BLOCK_LOG2;
        int b = _descriptor->_numBlocks;
        _allocator->allocate(&(_descriptor->_blocks[b]), nb - b);
        // zero out the newly allocated blocks
        for (; b < nb; ++b) {
            std::memset(
                map(_descriptor->_blocks[b] + (sizeof(ChunkEntryFlag) << ENTRIES_PER_BLOCK_LOG2)),
                0,
                sizeof(DataT) << ENTRIES_PER_BLOCK_LOG2
            );
        }
        _descriptor->_numBlocks = nb;
    }
}


/** Inserts the entry into the chunk, allocating memory if necessary. */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::insert(
    DataT          const & data,
    ChunkEntryFlag const   flags
) {
   int const   block = _descriptor->_nextBlock;
   std::size_t off   = _descriptor->_curBlockOffset;
   int         i     = _descriptor->_index;

   if (block == 0 || i >= (1 << ENTRIES_PER_BLOCK_LOG2)) {
       // no current block, or current block is full
       if (block >= _descriptor->_numBlocks) {
           if (block >= MAX_BLOCKS) {
               throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                   "maximum number of blocks per chunk exceeded");
           }
           off = _allocator->allocate();
           // zero out the newly allocated block
           std::memset(
               map(off + (sizeof(ChunkEntryFlag) << ENTRIES_PER_BLOCK_LOG2)),
               0,
               sizeof(DataT) << ENTRIES_PER_BLOCK_LOG2
           );
           _descriptor->_blocks[block] = off;
           _descriptor->_numBlocks = block + 1;
       } else {
           off = _descriptor->_blocks[block];
       }
       _descriptor->_nextBlock = block + 1;
       _descriptor->_curBlockOffset = off;
       i = 0;
   }

   // copy in data
   std::size_t const addr = off + i*sizeof(DataT) + 
                            (sizeof(ChunkEntryFlag) << ENTRIES_PER_BLOCK_LOG2);
   new (reinterpret_cast<DataT *>(map(addr))) DataT(data);
   *reinterpret_cast<ChunkEntryFlag *>(map(off + i*sizeof(ChunkEntryFlag))) = flags;
   _descriptor->_index = i + 1;
   ++_descriptor->_size;
}


/**
 * Walks through the chunk beginning at the @a i-th entry and removes all entries
 * marked as DELETED. Never throws.
 *
 * @param[in] i The index of the first entry to test for removal.
 * @return      @c true if any entries were removed
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
bool lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::pack(int const i) {

    static int const fb = (sizeof(ChunkEntryFlag) << ENTRIES_PER_BLOCK_LOG2);

    assert (i >= 0 && i < _descriptor->_size);

    int delta = 0x7fffffff;
    int sz    = _descriptor->_size;
    int dBlk  = i >> ENTRIES_PER_BLOCK_LOG2;
    int ib    = i & ((1 << ENTRIES_PER_BLOCK_LOG2) - 1);
    std::size_t d     = _descriptor->_blocks[dBlk];
    std::size_t dEnd  = d + BLOCK_SIZE;

    ChunkEntryFlag * df = map(d + ib*sizeof(ChunkEntryFlag));
    d += fb + ib*sizeof(DataT);

    int sBlk = dBlk;
    std::size_t s    = d;
    std::size_t sEnd = dEnd;
    ChunkEntryFlag const * sf = df;

    while (true) {
        for ( ; s < sEnd; s += sizeof(DataT), sf++) {
            if ((*sf & DELETED) != 0) {
                --sz;
                continue;
            }
            if ((*sf & IN_DELTA) != 0 && delta == 0x7fffffff) {
                // set delta
                delta = (dBlk << ENTRIES_PER_BLOCK_LOG2) +
                        (d - _descriptor->_blocks[dBlk] - fb)/sizeof(DataT);
            }
            if (df != sf) {
                *df = *sf;
                std::memcpy(map(d), map(s), sizeof(DataT));
            }
            df++;
            d += sizeof(DataT);
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
        std::size_t const cb = _descriptor->_blocks[dBlk];
        _descriptor->_curBlockOffset = cb;
        _descriptor->_nextBlock      = dBlk + 1;
        d -= cb + fb;
        _descriptor->_index = d/sizeof(DataT);
        _descriptor->_size  = sz;
        _descriptor->_delta = std::min(sz, delta);
        return true;
    }
    return false;
}


/**
 * Sets flag values for @a n entries in block @a b, starting at the @a i-th entry.
 *
 * @param[in] b     The block containing the entries for which flag values are to be set.
 * @param[in] flags The desired flag value of the block entries.
 * @param[in] i     The index of the first entry in block @a b for which the flag value is to be set.
 * @param[in] n     The number of block entries for which flag values are to be set.
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::setFlags(
    int const b,
    ChunkEntryFlag const flags,
    int const i,
    int const n
) {
    assert(i >= 0 && n >= 0 && i + n >= 0 && i + n <= (1 << ENTRIES_PER_BLOCK_LOG2));

    ChunkEntryFlag * df = getFlagBlock(b) + i;
    std::memset(df, flags, n);
}


/**
 * Walks through the chunk and undoes any uncommitted inserts or deletes. Never throws.
 *
 * @return  @c true if there were any modifications to rollback
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
bool lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::rollback() {
    bool mod = false;

    for (int b = 0; b < _descriptor->_nextBlock; ++b) {

        std::size_t const off = _descriptor->_blocks[b];
        ChunkEntryFlag * const flags = map(off);

        for (int i = 0, e = entries(b); i < e; ++i) {
            ChunkEntryFlag f = flags[i];
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


/**
 * Marks all uncommitted deletes/inserts as committed. Never throws.
 *
 * @param[in] clearDelta    If set to @c true, then the IN_DELTA flag bit is
 *                          cleared for each entry.
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::commit(bool clearDelta) {

    ChunkEntryFlag mask = UNCOMMITTED | INSERTED;
    if (clearDelta) {
        mask |= IN_DELTA;
    }
    mask = ~mask;

    for (int b = 0; b < _descriptor->_nextBlock; ++b) {
        std::size_t const off = _descriptor->_blocks[b];
        ChunkEntryFlag * const flags = map(off);
        int const e = entries(b);

        for (int i = 0; i < e; ++i) {
            flags[i] &= mask;
        }
    }

    if (clearDelta) {
        _descriptor->_delta = _descriptor->_size;
    }
}


/**
 * Marks the entries specified by the indexes in the given array as deleted.
 *
 * @param[in] deletes      An array containing the indexes of entries to delete.
 *                         Must be of length at least @a numDeletes.
 * @param[in] numDeletes   The number of entries to mark as deleted.
 * @param[in] end          Valid entry indexes must be less than this value.
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::applyDeletes(
    int const * const deletes,
    int const numDeletes,
    int const end
) {
    // check that all delete indexes are within the specified range
    for (int i = 0; i < numDeletes; ++i) {
        if (deletes[i] < 0 || deletes[i] >= end) {
            throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException,
                "Binary chunk delta file contains an invalid delete marker - delta not applied"
            );
        }
    }

    // ok - apply deletes.
    for (int i = 0; i < numDeletes; ++i) {
        int d = deletes[i];
        ChunkEntryFlag * f = reinterpret_cast<ChunkEntryFlag *>(map(
            _descriptor->_blocks[d >> ENTRIES_PER_BLOCK_LOG2] +
            (d & ((1 << ENTRIES_PER_BLOCK_LOG2) - 1)) * sizeof(ChunkEntryFlag)
        ));
        *f |= DELETED;
    }
}


namespace lsst { namespace ap { namespace {

void doRead(io::SequentialReader & reader, unsigned char * dst, std::size_t dstlen) {
     while (dstlen > 0) {
        std::size_t nb = reader.read(dst, dstlen);
        assert(nb <= dstlen);
        if (nb == 0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException, "Unexpected end of file");
        }
        dst    += nb;
        dstlen -= nb;
    }
}

std::string const badChunkMessage =
    "Binary chunk file failed sanity checks - file is of the wrong format, "
    "of the wrong record type, or corrupted. This error indicates that one "
    "of the following situations has occurred:\n"
    "  - the file being read was not a chunk file\n"
    "  - the chunk file being read was generated with an incompatible\n"
    "    version of the association pipeline code\n"
    "  - the chunk file being read was generated on a machine with\n"
    "    either a different word size (32 vs. 64 bit) or endianness from\n"
    "    the runtime machine(s).\n"
    "Please regenerate the chunk on a machine matching the characteristics "
    "of the runtime machine(s), check your policy files, and try again.";

}}} // end of namespace lsst::ap::<anonymous>


/**
 * Reads the data from the binary chunk file @a name into this chunk. Note that this chunk is
 * emptied immediately on entering the function.
 *
 * @param name         The name of binary chunk file to read into memory
 * @param compressed   Is the binary chunk file compressed? zlib or gzip compression is supported.
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::read(
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
    doRead(*reader, reinterpret_cast<unsigned char *>(&header), sizeof(BinChunkHeader));
    if (!header.isValid() || header._numDeletes != 0 || header._recordSize != sizeof(DataT)) {
        throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException, badChunkMessage);
    }

    int nr = header._numRecords;
    if (nr == 0) {
        return; // nothing to read in
    }
    reserve(nr);

    int b = 0;
    int nd;

    // read in one memory block at a time
    do {
        nd  = std::min((1 << ENTRIES_PER_BLOCK_LOG2), nr);
        nr -= nd;
        doRead(*reader, reinterpret_cast<unsigned char *>(getBlock(b)), nd * sizeof(DataT));
        setFlags(b, 0, 0, nd);
        ++b;
    } while (nr > 0);

    // update chunk size
    _descriptor->_nextBlock = b;
    _descriptor->_curBlockOffset = _descriptor->_blocks[b - 1];
    _descriptor->_index = nd;
    int sz = nd + ((b - 1) << ENTRIES_PER_BLOCK_LOG2);
    _descriptor->_size  = sz;
    _descriptor->_delta = sz;
}


/**
 * Reads the given binary chunk delta file, performing any indicated deletes and appending new
 * records to the end of this chunk. This function provides the strong exception safety
 * guarantee.
 *
 * @param name         The name of binary chunk file to read into memory
 * @param compressed   Is the binary chunk file compressed? zlib or gzip compression is supported.
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::readDelta(
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
    doRead(*reader, reinterpret_cast<unsigned char *>(&header), sizeof(BinChunkHeader));
    if (!header.isValid() || header._recordSize != sizeof(DataT)) {
        throw LSST_EXCEPT(lsst::pex::exceptions::IoErrorException, badChunkMessage);
    }

    // read in indexes of records to delete
    boost::scoped_array<int> deletes(header._numDeletes > 0 ? new int[header._numDeletes] : 0);
    doRead(*reader, reinterpret_cast<unsigned char *>(deletes.get()), sizeof(int) * header._numDeletes);

    // read in records to append
    int nr = header._numRecords;
    if (nr == 0) {
        applyDeletes(deletes.get(), header._numDeletes, _descriptor->_size);
        return; // nothing more to read in
    }
    reserve(nr + _descriptor->_size);

    // fill up the current block (or the first block if there is no current block)
    int b  = _descriptor->_nextBlock;
    if (b > 0) {
        --b;
    }
    int i = _descriptor->_index;
    int nd = std::min((1 << ENTRIES_PER_BLOCK_LOG2) - i, nr);
    nr -= nd;
    doRead(*reader, reinterpret_cast<unsigned char *>(&getBlock(b)[i]), nd * sizeof(DataT));
    setFlags(b, IN_DELTA, i, nd);
    nd += i;
    ++b;

    // fill remaining blocks
    while (nr > 0) {
        nd  = std::min((1 << ENTRIES_PER_BLOCK_LOG2), nr);
        nr -= nd;
        doRead(*reader, reinterpret_cast<unsigned char *>(getBlock(b)), nd * sizeof(DataT));
        setFlags(b, IN_DELTA, 0, nd);
        ++b;
    }

    int const sz = ((b - 1) << ENTRIES_PER_BLOCK_LOG2) + nd;
    applyDeletes(deletes.get(), header._numDeletes, sz);

    // update chunk state to reflect additions
    _descriptor->_nextBlock = b;
    _descriptor->_curBlockOffset = _descriptor->_blocks[b - 1];
    _descriptor->_index = nd;
    _descriptor->_size = sz;
}


/**
 * Writes the data from this chunk to a binary file. Note that deleted records will be
 * written out as well - to skip them, call pack() immediately before this function. Note
 * that uncomitted deletes/inserts and entries flagged as IN_DELTA do not automatically
 * have their status cleared - to do this call commit() with @c true as the argument.
 *
 * @param name         The name of binary chunk file to write.
 * @param overwrite    Should an existing file with the given name be overwritten?
 * @param compressed   Should the binary chunk file contents be compressed? A gzip compatible
 *                     algorithm will be used.
 * @param withDelta    Should entries marked IN_DELTA be written out?
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::write(
    std::string const & name,
    bool        const   overwrite,
    bool        const   compressed,
    bool        const   withDelta
) const {
    boost::scoped_ptr<io::SequentialWriter> writer;
    if (compressed) {
        writer.reset(new io::CompressedFileWriter(name, overwrite));
    } else {
        writer.reset(new io::SequentialFileWriter(name, overwrite));
    }

    // create header
    int nr = withDelta ? size() : delta();
    BinChunkHeader header;
    header._numRecords = nr;
    header._recordSize = sizeof(DataT);

    writer->write(reinterpret_cast<unsigned char *>(&header), sizeof(BinChunkHeader));
    for (int b = 0; nr > 0; ++b) {
        int const nd = std::min((1 << ENTRIES_PER_BLOCK_LOG2), nr);
        writer->write(reinterpret_cast<unsigned char const *>(getBlock(b)), nd * sizeof(DataT));
        nr -= nd;
    }
    writer->finish();
}


/**
 * Writes any deletes and inserts in this chunk to a binary delta file named @a name. Note that
 * even on successful function return, uncommitted deletes/inserts are @b not marked as committed
 * (allowing for writeDelta calls on multiple chunks to be aggregated into a single transaction).
 * Finally, note that deleted delta records are written out by this function. To skip them, call
 * pack() on this chunk, passing delta() as the parameter.
 *
 * @param name         The name of binary chunk file to write
 * @param overwrite    Should an existing file with the given name be overwritten?
 * @param compressed   Should the binary chunk delta file contents be compressed (a gzip
 *                     compatible algorithm will be used)?
 */
template <typename AllocatorT, typename DataT, typename TraitsT>
void lsst::ap::ChunkRef<AllocatorT, DataT, TraitsT>::writeDelta(
    std::string const & name,
    bool        const   overwrite,
    bool        const   compressed
) const {
    boost::scoped_ptr<io::SequentialWriter> writer;
    if (compressed) {
        writer.reset(new io::CompressedFileWriter(name, overwrite));
    } else {
        writer.reset(new io::SequentialFileWriter(name, overwrite));
    }

    // collect all deletes (there are likely to be very few, if any)
    std::vector<int> deletes;

    int const nb = _descriptor->_nextBlock;
    for (int b = 0; b < nb; ++b) {
        ChunkEntryFlag const * const f = getFlagBlock(b);
        int const nd = (b < nb - 1) ? (1 << ENTRIES_PER_BLOCK_LOG2) : _descriptor->_index;
        for (int i = 0; i < nd; ++i) {
            if ((f[i] & DELETED) != 0) {
                deletes.push_back(i + (b << ENTRIES_PER_BLOCK_LOG2));
            }
        }
    }

    // create and write header
    int fd = _descriptor->_delta;
    BinChunkHeader header;
    header._numRecords = _descriptor->_size - fd;
    header._numDeletes = deletes.size();
    header._recordSize = sizeof(DataT);
    writer->write(reinterpret_cast<unsigned char *>(&header), sizeof(BinChunkHeader));

    // write array of delete indexes
    if (!deletes.empty()) {
        writer->write(reinterpret_cast<unsigned char *>(&deletes.front()),
                      deletes.size() * sizeof(int));
    }

    // write out IN_DELTA chunk entries
    if (fd < _descriptor->_size) {
        int b = fd >> ENTRIES_PER_BLOCK_LOG2;
        fd &= (1 << ENTRIES_PER_BLOCK_LOG2) - 1;
        int nd = (b < nb - 1) ? (1 << ENTRIES_PER_BLOCK_LOG2) : _descriptor->_index;
        writer->write(reinterpret_cast<unsigned char const *>(&getBlock(b)[fd]), (nd - fd) * sizeof(DataT));
        for (++b; b < nb; ++b) {
            nd = (b < nb - 1) ? (1 << ENTRIES_PER_BLOCK_LOG2) : _descriptor->_index;
            writer->write(reinterpret_cast<unsigned char const *>(getBlock(b)), nd * sizeof(DataT));
        }
    }

    // all done
    writer->finish();
}


#endif // LSST_AP_CHUNK_CC
