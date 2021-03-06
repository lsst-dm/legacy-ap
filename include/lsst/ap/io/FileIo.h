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
 * @brief   Low-level sequential file IO classes.
 *
 * @ingroup ap
 */

#ifndef LSST_AP_IO_FILE_IO_H
#define LSST_AP_IO_FILE_IO_H

#include <aio.h>
#include <zlib.h>

#include <string>

#include "boost/noncopyable.hpp"
#include "boost/scoped_array.hpp"

#include "../Common.h"


namespace lsst { namespace ap { namespace io {

/** @brief  Abstract base class for sequential I/O classes. */
class SequentialIoBase {

public :

    enum State {
        IN_PROGRESS = 0,
        FINISHED,
        FAILED
    };

    SequentialIoBase();

    virtual ~SequentialIoBase() = 0;

    /// Returns @c true if there are no more bytes available for reading.
    bool finished() const { return _state == FINISHED; }

    /// Returns @c true if a read operation failed.
    bool failed() const { return _state == FAILED; }

    /// Returns the state of the SequentialReader.
    State getState() const { return _state; }

protected :

    State _state;
};


/** @brief  Abstract base class for reading a stream of data in sequential fashion. */
class SequentialReader : public SequentialIoBase {

public :

    /**
     * Reads up to @a len bytes from an underlying storage device into @a buf
     * and returns the number of bytes actually read.
     */
    virtual std::size_t read(unsigned char * const buf, std::size_t const len) = 0;
};


/** @brief  Abstract base class for writing a stream of data in sequential fashion. */
class SequentialWriter : public SequentialIoBase {

public :

    /// Writes @a len bytes from @a buf to the underlying storage device.
    virtual void write(unsigned char const * const buf, std::size_t const len) = 0;

    /// Moves modified data to the underlying storage device and marks the SequentialWriter as finished.
    virtual void finish() = 0;
};


/** @brief  A sequential reader for uncompressed files. Uses standard (blocking) IO calls. */
class SequentialFileReader :
    public  SequentialReader,
    private boost::noncopyable
{

public:

    explicit SequentialFileReader(std::string const & fileName);

    virtual ~SequentialFileReader();
    virtual std::size_t read(unsigned char * const buf, std::size_t const size);

private :

    int _fd;

    void cleanup();
    void cleanup(State const state) {
        cleanup();
        _state = state;
    }
};


/** @brief  A sequential writer for uncompressed files. Uses standard (blocking) IO calls. */
class SequentialFileWriter :
    public  SequentialWriter,
    private boost::noncopyable
{

public :

    explicit SequentialFileWriter(
        std::string const & fileName,
        bool        const   overwrite = false
    );

    virtual ~SequentialFileWriter();
    virtual void write(unsigned char const * const buf, std::size_t const len);
    virtual void finish();

private :

    int _fd;

    void cleanup();
    void cleanup(State const state) {
        cleanup();
        _state = state;
    }
};


/**
 * @brief   A sequential reader for compressed files that uses asynchronous IO
 *          to overlap IO with decompression.
 *
 * Gzip compatible files, or files written with zlib compression can be read by this class.
 */
class CompressedFileReader : public SequentialReader {

public :

    explicit CompressedFileReader(
        std::string const & fileName,
        std::size_t const   blockSize = 262144
    );

    virtual ~CompressedFileReader();
    virtual std::size_t read(unsigned char * const buf, std::size_t const size);

    std::size_t getBlockSize() const { return _blockSize; }

private :

    boost::scoped_array<unsigned char> _memory;
    unsigned char *   _buffers;   ///< aligned input buffers
    ::z_stream        _stream;    ///< zlib state
    ::aiocb           _request;   ///< Outstanding IO request
    std::size_t const _blockSize; ///< read granularity
    std::size_t       _fileSize;  ///< Size of the file being read
    std::size_t       _remaining; ///< Bytes that haven't yet been read
    int               _fd;        ///< file descriptor

    void cleanup();
    void cleanup(State const state) {
        cleanup();
        _state = state;
    }
};


/**
 * @brief   A sequential writer for compressed files that uses asynchronous IO
 *          to overlap IO with compression.
 *
 * Gzip compatible files are written by this class.
 */
class CompressedFileWriter : public SequentialWriter {

public :

    explicit CompressedFileWriter(
        std::string const & fileName,
        bool        const   overwrite = false,
        std::size_t const   blockSize = 262144
    );

    virtual ~CompressedFileWriter();

    virtual void write(unsigned char const * const buf, std::size_t const size);
    virtual void finish();

    std::size_t getBlockSize() const { return _blockSize; }

private :

    boost::scoped_array<unsigned char> _memory;
    unsigned char *   _buffers;   ///< aligned output buffers
    ::z_stream        _stream;    ///< zlib state
    ::aiocb           _request;   ///< Outstanding IO request
    std::size_t const _blockSize; ///< read granularity
    int               _fd;
    bool              _started;

    void cleanup();
    void cleanup(State const state) {
        cleanup();
        _state = state;
    }
};


}}}  // end of namespace lsst::ap::io

#endif // LSST_AP_IO_FILE_IO_H
